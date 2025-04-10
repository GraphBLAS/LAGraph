//------------------------------------------------------------------------------
// LAGraph_RichClubCoefficient: rich club coefficient of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel Gomez, Texas A&M University

//------------------------------------------------------------------------------

// Get the rich club coefficient of a graph.

// Given a Symetric Graph with no self edges, LAGraph_RichClubCoefficient will
// calculate the rich club coefficients of the graph. 

// The values will be output as a sparse GrB_Vector, the rich club coefficient 
// of k will be found at the closest entry at or above k.

// The G->out_degree cached property must be defined for this method.

// References:

// Julian J. McAuley, Luciano da Fontoura Costa, and Tibério S. Caetano, “The 
// rich-club phenomenon across complex network hierarchies”, Applied Physics 
// Letters Vol 91 Issue 8, August 2007. https://arxiv.org/abs/physics/0701290

#define LG_FREE_WORK                                    \
{                                                       \
    /* free any workspace used here */                  \
    GrB_free(&D) ;                                      \
    GrB_free(&P) ;                                      \
    GrB_free(&A_deg) ;                                  \
    GrB_free(&degrees) ;                                \
    GrB_free(&deg_x) ;                                  \
    GrB_free(&node_edges) ;                             \
    GrB_free(&node_edges_x) ;                           \
    GrB_free(&ones_v) ;                                 \
    GrB_free(&edges_per_deg) ;                          \
    GrB_free(&verts_per_deg) ;                          \
    GrB_free(&iseq_2lt) ;                               \
    GrB_free(&plus_2le) ;                               \
    GrB_free(&rcCalculation) ;                          \
    GrB_free(&ramp_v) ;                                 \
    LAGraph_Free(&a_space, NULL) ;                      \
}


#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (rccs) ;      \
}

#include "LG_internal.h"
#include "LAGraphX.h"

typedef void (*LAGraph_binary_function) (void *, const void *, const void *) ;

#define ISEQ_2ISLT                                                          \
    "void iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y)            \n"\
    "{                                                                          \n"\
        "(*z) = (int64_t)((*x < *y) + (*x <= *y)) ;                             \n"\
    "}"
void iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y)
{
    (*z) = (int64_t)((*x < *y) + (*x <= *y)) ;
}

#define RICH_CLUB_FORMULA                                                      \
    "void rich_club_formula(double *z, const int64_t *x, const int64_t *y)      \n"\
    "{                                                                          \n"\
    "   (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0)) ;     \n"\
    "}"
void rich_club_formula(double *z, const int64_t *x, const int64_t *y)
{
    (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0));
} 
int LAGraph_RichClubCoefficient
(
    // output:
    //rccs(i): rich club coefficent of i
    GrB_Vector *rccs,    

    // input: 
    LAGraph_Graph G, //input graph
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    LG_CLEAR_MSG ;

    // n x n Adjacency Matrix
    // With values cooresponding to the degree of its column
    GrB_Matrix A_deg = NULL;

    // n x n Diagonal Matrix
    // entries corresponding to degrees.
    GrB_Matrix D = NULL;

    // n degrees vector
    GrB_Vector degrees = NULL, deg_x = NULL;

    // n x 1
    // contains the number of edges for which the ith node is
    // the smallest degree node * 2 + # edges w/ same degree as the other node
    // to account for double counting of edges w/ same degree as the other node.
    GrB_Vector node_edges = NULL, node_edges_x = NULL;

    // max_degree x 1
    // the ith entry contains the number of edges whose lowest degree is i.
    GrB_Vector edges_per_deg = NULL;

    // max_degree x 1
    // the ith entry contains the number of verticies whose degree is i.
    GrB_Vector verts_per_deg = NULL;

    // edge_vec_nvals x 1
    // Vector of ones
    GrB_Vector ones_v = NULL;

    // Ramp vector
    GrB_Vector ramp_v = NULL;

    // 2 * (x < y) + (x == y)
    GrB_BinaryOp iseq_2lt = NULL;

    // [+].[iseq_2lt]
    GrB_Semiring plus_2le = NULL;

    // 2E_K / (N_k (N_k -1))
    GrB_BinaryOp rcCalculation = NULL;

    GrB_Matrix A = NULL; // G->A, the adjacency matrix

    // Matrix used for row reduction
    GrB_Matrix P = NULL;

    GrB_Index n ;
    
    GrB_Index edge_vec_nvals;
    GrB_Index max_deg;
    bool iso = false;

    void *a_space = NULL;
    
    int64_t *node_edges_arr = NULL, *deg_arr = NULL, 
        *epd_arr = NULL, *ones = NULL, 
        *vpd_arr = NULL;
    GrB_Type epd_type = NULL, vpd_type = NULL;
    int64_t epd_n = 0, vpd_n = 0, epd_size = 0, vpd_size = 0;
    int epd_h = 0, vpd_h = 0;
    GrB_Index *epd_index = NULL,  *vpd_index = NULL;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rccs != NULL, GrB_NULL_POINTER);

    LG_ASSERT_MSG(
        G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE, 
        "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,
        "G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE, 
        "G->nself_edges must be zero") ; 

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new(&A_deg, GrB_INT64,n,n)) ;

    GRB_TRY (GrB_Vector_new(&degrees, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_new(&node_edges, GrB_INT64, n)) ;
    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_BinaryOp_new(
        &iseq_2lt, (LAGraph_binary_function) (&iseq_2islt), 
        GrB_INT64, GrB_INT64, GrB_INT64, "iseq_2islt", ISEQ_2ISLT)) ;
    GRB_TRY (GxB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&rich_club_formula), 
        GrB_FP64, GrB_INT64, GrB_INT64, 
        "rich_club_formula", RICH_CLUB_FORMULA)) ;
    #else
    GRB_TRY (GrB_BinaryOp_new(
        &iseq_2lt, (LAGraph_binary_function) (&iseq_2islt), 
        GrB_INT64, GrB_INT64, GrB_INT64)) ;
    GRB_TRY (GrB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&rich_club_formula), 
        GrB_FP64, GrB_INT64, GrB_INT64 )) ;
    #endif

    GRB_TRY (GrB_Semiring_new(&plus_2le, GrB_PLUS_MONOID_INT64, iseq_2lt)) ;
    
    GRB_TRY (GrB_Vector_reduce_INT64(
        &max_deg, NULL, GrB_MAX_MONOID_INT64, G->out_degree, NULL)) ;
    GRB_TRY (GrB_Vector_new(&edges_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY (GrB_Vector_new(&verts_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY (GrB_Vector_new(rccs, GrB_FP64, max_deg)) ;

    //--------------------------------------------------------------------------
    // Calculations
    //--------------------------------------------------------------------------

    // degrees = G->out_degree - 1
    // Fill out degree vector, to target col_scale mxm on graphs 
    // with singletons, scalar value irrelevant.
    GRB_TRY (GrB_Vector_assign_INT64(
        degrees, NULL, NULL, (int64_t) -1, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign(
        degrees, NULL, GrB_PLUS_INT64, G->out_degree, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_diag(&D, degrees, 0)) ;

    // Each edge in the graph gets the value of the degree of its row node
    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GrB_mxm(
        A_deg, NULL, NULL, GxB_ANY_FIRST_INT64, D, A, NULL)) ;
    #else
    GRB_TRY (GrB_mxm(
        A_deg, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, D, A, NULL)) ;
    #endif
    // Sum the number of edges each node is "responsible" for.
    GRB_TRY (GrB_mxv(
        node_edges, NULL, GrB_PLUS_INT64, plus_2le, A_deg, degrees, NULL)) ;

    // The rest of this is indexing the number of edges and number of nodes at 
    // each degree and then doing a cummulative sum to know the amount of edges 
    // and nodes at degree geq k.
    GRB_TRY (GrB_Vector_nvals (&edge_vec_nvals, node_edges)) ;
    #if USING_GRAPHBLAS_V10
        if(n == edge_vec_nvals)
        {
            deg_x = degrees;
            degrees = NULL;
            node_edges_x = node_edges;
            node_edges = NULL;
        }
        else
        {
            GRB_TRY (GrB_Vector_assign(
                degrees, G->out_degree, NULL, degrees, GrB_ALL, 0, GrB_DESC_RS
            )) ;
            GRB_TRY (GrB_Vector_new(&deg_x, GrB_BOOL, 0)) ;  
            GRB_TRY (GrB_Vector_new(&node_edges_x, GrB_BOOL, 0)) ;  
            GRB_TRY (GxB_Vector_extractTuples_Vector(
                NULL, deg_x, degrees, NULL
            )) ;
            GRB_TRY (GxB_Vector_extractTuples_Vector(
                NULL, node_edges_x, node_edges, NULL
            )) ;
        }
        GRB_TRY (GrB_Vector_nvals(&edge_vec_nvals, node_edges_x))
        GRB_TRY (GrB_Vector_new(&ones_v, GrB_INT64, edge_vec_nvals)) ;
        GRB_TRY (GrB_Vector_new(&ramp_v, GrB_INT64, edge_vec_nvals + 1)) ;  

        GRB_TRY (GrB_Vector_assign_INT64(
            ramp_v, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            edges_per_deg, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            verts_per_deg, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            ones_v, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;

        GRB_TRY (GrB_apply (
            ramp_v, NULL, NULL, GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
        LG_TRY (LAGraph_FastAssign (
            edges_per_deg, NULL, GrB_PLUS_INT64, deg_x, node_edges_x, ramp_v,
            GxB_PLUS_SECOND_INT64, NULL, msg
        )) ;
        LG_TRY (LAGraph_FastAssign (
            verts_per_deg, NULL, GrB_PLUS_INT64, deg_x, ones_v, ramp_v,
            GxB_PLUS_PAIR_INT64, NULL, msg
        )) ;

        GRB_TRY (GxB_Vector_unload(
            edges_per_deg, (void **) &epd_arr, &epd_type,
            &epd_n, &epd_size, &epd_h, NULL)) ;
        GRB_TRY (GxB_Vector_unload(
            verts_per_deg, (void **) &vpd_arr, &vpd_type,
            &vpd_n, &vpd_size, &vpd_h, NULL)) ;
        
        LG_ASSERT (max_deg == vpd_n && max_deg == epd_n, GrB_INVALID_VALUE) ;
        //run a cummulative sum (backwards) on vpd_arr
        for(GrB_Index i = max_deg - 1; i > 0; --i)
        {
            vpd_arr[i-1] += vpd_arr[i] ;
            epd_arr[i-1] += epd_arr[i] ;
        }
        GRB_TRY(GxB_Vector_load(
            edges_per_deg, (void **) &epd_arr, epd_type,
            epd_n, epd_size, epd_h, NULL)) ;
        GRB_TRY(GxB_Vector_load(
            verts_per_deg, (void **) &vpd_arr, vpd_type,
            vpd_n, vpd_size, vpd_h, NULL)) ;
    #else
        LG_TRY (LAGraph_Malloc(
            &a_space, edge_vec_nvals * 3 + max_deg * 4, sizeof(int64_t), NULL
        )) ;
        int64_t *T = a_space;
        deg_arr = T;            T += edge_vec_nvals;
        node_edges_arr = T;     T += edge_vec_nvals;
        ones = T;               T += edge_vec_nvals;
        epd_arr = T;            T += max_deg;
        vpd_arr = T;            T += max_deg;
        epd_index = T;          T += max_deg;
        vpd_index = T;          T += max_deg;

        #pragma omp parallel for schedule(static)
        for(uint64_t i = 0; i < edge_vec_nvals; ++i)
        {
            ones[i] = 1ll;
        }
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64(
            degrees, NULL, NULL, GrB_MINUS_INT64, G->out_degree, 1, NULL)) ;
        //TODO: remove NULL for Vanilla GB
        GRB_TRY (GrB_Vector_extractTuples_INT64(
            NULL, deg_arr, &edge_vec_nvals, degrees
        )) ;
        GRB_TRY (GrB_Vector_extractTuples_INT64(
            NULL, node_edges_arr, &edge_vec_nvals, node_edges
        )) ;

        // Build with degrees as indecies and handle duplicates via adition
        GRB_TRY (GrB_Vector_build_INT64 (
            edges_per_deg, deg_arr, node_edges_arr, edge_vec_nvals, 
            GrB_PLUS_INT64)) ;
        GRB_TRY (GrB_Vector_build_INT64 (
            verts_per_deg, deg_arr, ones, edge_vec_nvals, GrB_PLUS_INT64)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            edges_per_deg, edges_per_deg, NULL, (int64_t) 0, 
            GrB_ALL, 0, GrB_DESC_SC)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            verts_per_deg, verts_per_deg, NULL, (int64_t) 0, 
            GrB_ALL, 0, GrB_DESC_SC)) ;
        
        // Extract into arrays
        GRB_TRY (GrB_Vector_extractTuples_INT64(
            epd_index, epd_arr, &max_deg, edges_per_deg
        )) ;
        GRB_TRY (GrB_Vector_extractTuples_INT64(
            vpd_index, vpd_arr, &max_deg, verts_per_deg
        )) ;
        //run a cummulative sum (backwards) on vpd_arr
        for(GrB_Index i = max_deg - 1; i > 0; --i)
        {
            vpd_arr[i-1] += vpd_arr[i] ;
            epd_arr[i-1] += epd_arr[i] ;
        }
        GRB_TRY (GrB_Vector_clear(edges_per_deg)) ;
        GRB_TRY (GrB_Vector_clear(verts_per_deg)) ;
        GRB_TRY (GrB_Vector_build_INT64(
            edges_per_deg, epd_index, epd_arr, max_deg, NULL
        )) ;
        GRB_TRY (GrB_Vector_build_INT64(
            verts_per_deg, vpd_index, vpd_arr, max_deg, NULL
        )) ;
        T = deg_arr = node_edges_arr = ones = NULL ;
        epd_index = vpd_index = epd_arr = vpd_arr = NULL ;
    #endif

    /**
     * Cumulative sum (TODO: should be a GBLAS method!)
     * 
     * GrB_cumsum(GrB_Matrix C, const GrB_Matrix mask, const GrB_BinaryOp accum,
     *      const GrB_BinaryOp plus, GrB_Matrix A, const GrB_Descriptor desc)
     * 
     * By default sums rows. Returns a nearly full matrix:
     * [., ., 1, 1, 1, 1, ., ., 1] --> [., ., 1, 2, 3, 4, 4, 4, 5]
     * Mask can be A, then returns a matrix with the same pattern.
     * [., ., 1, 1, 1, 1, ., ., 1] --> [., ., 1, 2, 3, 4, ., ., 5]
     * 
     * Should we be able to sum in the opposite direction? 
     *  Yes since not all monoids have inverse operations. 
     * 
     * If plus biop is not a monoid, this method should still work?
     */
    
    //Computes the RCC of a matrix
    GRB_TRY(GrB_eWiseMult(
        *rccs, NULL, NULL, rcCalculation, 
        edges_per_deg, verts_per_deg, NULL
    )) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
#undef LG_FREE_WORK
#undef LG_FREE_ALL
#define LG_FREE_WORK                                    \
{                                                       \
    /* free any workspace used here */                  \
    LAGraph_Free(&a_space, NULL) ;                  \
    GrB_free (&cont) ;                                  \
    LAGraph_Free((void **)&Ai, NULL) ;               \
    LAGraph_Free((void **)&Ap, NULL) ;               \
    LAGraph_Free((void **)&slice, NULL) ;               \
}


#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    LAGraph_Free((void **)&rcc, NULL) ;     \
    GrB_free (rccs) ;      \
}
#define TIMINGS
#ifdef TIMINGS
static void print_timings (const double timings [16])
{
    double total = timings [0] + timings [1] + timings [2] + timings [3] + timings [4];
    printf ("RCC %12.6f (%4.1f%%) init\n", timings [0], 100. * timings [0] / total) ;
    printf ("RCC %12.6f (%4.1f%%) counting edges\n", timings [1], 100. * timings [1] / total) ;
    printf ("RCC %12.6f (%4.1f%%) counting nodes\n", timings [2], 100. * timings [2] / total) ;
    printf ("RCC %12.6f (%4.1f%%) cumulative sum\n", timings [3], 100. * timings [3] / total) ;
    printf ("RCC %12.6f (%4.1f%%) calculation\n", timings [4], 100. * timings [4] / total) ;
}
#endif

//Scuffed upperbound function 
static int64_t LG_binary_search    // returns upperbound - 1
(
    const int64_t pivot,
    const int64_t *LG_RESTRICT X_0,         // search in X [p_start..p_end_-1]
    const int64_t p_start,
    const int64_t p_end
)
{

    //--------------------------------------------------------------------------
    // find where the Pivot appears in X
    //--------------------------------------------------------------------------

    // binary search of X [p_start...p_end-1] for the Pivot
    int64_t pleft = p_start ;
    int64_t pright = p_end;
    while (pleft < pright)
    {
        int64_t pmiddle = pleft + (pright - pleft) / 2 ;
        bool less = (X_0 [pmiddle] < pivot) ;
        pleft  = less ? pmiddle + 1 : pleft ;
        pright = less ? pright : pmiddle ;
    }
    if(X_0[pleft] <= pivot)
        pleft++;
    return (--pleft) ;
}


int LAGraph_RichClubCoefficient_NoGB
(
    // output:
    //rccs(i): rich club coefficent of i
    GrB_Vector *rccs,    

    // input: 
    LAGraph_Graph G, //input graph
    char *msg
)
{
    #if USING_GRAPHBLAS_V10
    GxB_Container cont = NULL;
    GrB_Matrix A = G->A;
    int64_t  *Ap = NULL, *Ai = NULL;
    void *a_space = NULL;
    GrB_Type p_type = NULL, i_type = NULL;
    int p_hand = 0, i_hand = 0;
    int n_threads = LG_nthreads_outer * LG_nthreads_inner;
    uint64_t p_n = 0, i_n = 0, p_size = 0, i_size = 0, max_deg = 0;
    uint64_t *epd = NULL, *vpd = NULL;
    int64_t *LG_RESTRICT slice  = NULL;
    double *rcc = NULL;
    #ifdef TIMINGS
    double timings [16] ;
    memset(timings, 0, 16*sizeof(double)) ;
    double tic = LAGraph_WallClockTime ( ) ;
    LG_SET_BURBLE (false) ;
    #endif

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rccs != NULL, GrB_NULL_POINTER);

    LG_ASSERT_MSG(
        G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE, 
        "G->A must be symmetric") ;
    LG_ASSERT_MSG(
        G->is_symmetric_structure == LAGraph_TRUE, GrB_INVALID_VALUE, 
        "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,
        "G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE, 
        "G->nself_edges must be zero") ; 
    GRB_TRY(GxB_Container_new(&cont)) ;
    GRB_TRY(GxB_unload_Matrix_into_Container(A, cont, NULL)) ;
    LG_ASSERT_MSG(cont->format == GxB_SPARSE, GrB_NOT_IMPLEMENTED, 
        "Matrix must be sparse") ;    
    LG_TRY (LAGraph_Malloc(
        (void **) &Ap, cont->nvals, sizeof(uint64_t), NULL)) ;
    LG_TRY (LAGraph_Malloc(
        (void **) &Ai, cont->nvals + 1, sizeof(uint64_t), NULL)) ;
    p_n = cont->nvals + 1; i_n = cont->nvals;
    GRB_TRY (GrB_Vector_extractTuples_INT64(
        NULL, Ap, &p_n, cont->p)) ;
    GRB_TRY (GrB_Vector_extractTuples_INT64(
        NULL, Ai, &i_n, cont->i)) ;
    GRB_TRY (GxB_load_Matrix_from_Container(A, cont, NULL)) ;
    GRB_TRY (GrB_Vector_reduce_INT64(
        &max_deg, NULL, GrB_MAX_MONOID_INT64, G->out_degree, NULL)) ;
    int64_t i = 0;
    #ifdef TIMINGS
    timings[0] = LAGraph_WallClockTime ( );
    #endif
    LG_TRY (LAGraph_Calloc(&a_space, max_deg * 2, sizeof(uint64_t), NULL)) ;
    LG_TRY (LAGraph_Malloc((void **)&slice, n_threads + 1, sizeof(int64_t), NULL)) ;
    epd = a_space ;
    vpd = a_space + max_deg * sizeof(uint64_t) ;
    LG_TRY (LAGraph_Malloc((void **) &rcc, max_deg, sizeof(double), NULL)) ;
    // while(ptr < p_n - 1)
    // {
    //     uint64_t dp = Ap[ptr+1] - Ap[ptr];
    //     for(; Ap[ptr + 1] > i; ++i)
    //     {
    //         uint64_t di = Ap[Ai[i]+1] - Ap[Ai[i]] ;
    //         epd[dp - 1] += (dp < di) + (dp <= di) ;
    //     }
    //     if (dp > 0)
    //         ++vpd[dp - 1] ;
    //     ++ptr ;
    // }
    LG_eslice (slice, i_n, n_threads) ;
    #pragma omp parallel for num_threads(n_threads) schedule(static, 1) private(i)
    for (int tid = 0 ; tid < n_threads ; tid++)
    {
        int64_t loc_sum = 0, dp = 0;
        int64_t loc_arr[1024];
        memset(loc_arr, 0, 1024 * sizeof(int64_t));
        i = slice[tid];
        int64_t ptr = LG_binary_search(i, Ap, 0, p_n - 1) ;
        while(i < slice[tid + 1])
        {
            while(Ap[ptr + 1] <= i) ++ptr;
            int64_t dp = Ap[ptr + 1] - Ap[ptr];
            if(dp <= 1024)
                for(; i < slice[tid + 1] && i < Ap[ptr + 1]; ++i)
                {
                    uint64_t di = Ap[Ai[i] + 1] - Ap[Ai[i]];
                    loc_arr[dp - 1] += (dp < di) + (dp <= di);
                }
            else
            {
                loc_sum = 0;
                for(; i < slice[tid + 1] && i < Ap[ptr + 1]; ++i)
                {
                    uint64_t di = Ap[Ai[i]+1] - Ap[Ai[i]];
                    loc_sum += (dp < di) + (dp <= di);
                }
                #pragma omp atomic
                    epd[dp - 1] += loc_sum ;
            }
        }
        #pragma omp critical
        {
            for(int64_t j = 0; j < 1024 && j < max_deg; ++j)
            {
                epd[j] += loc_arr[j];
            }
        }
    }
    #ifdef TIMINGS
    timings[1] = LAGraph_WallClockTime ( );
    #endif
    
    #pragma omp parallel
    {
        int64_t loc_arr[1024];
        memset(loc_arr, 0, 1024 * sizeof(int64_t));
        #pragma omp for schedule(static)
        for(i = 0; i < p_n - 1; ++i)
        {
            int64_t dp = Ap[i + 1] - Ap[i] - 1;
            if(dp < 0) continue;
            if(dp < 1024)
            {
                ++loc_arr[dp];
            }
            else
            {
                #pragma omp atomic
                    ++vpd[dp];
            }
        }  
        #pragma omp critical
        {
            for(int64_t j = 0; j < 1024 && j < max_deg; ++j)
            {
                vpd[j] += loc_arr[j];
            }
        }
    }
    
    #ifdef TIMINGS
    timings[2] = LAGraph_WallClockTime ( );
    #endif
    //run a cummulative sum (backwards)
    for(i = max_deg - 1; i > 0; --i)
    {
        vpd[i-1] += vpd[i] ;
        epd[i-1] += epd[i] ;
    }
    #ifdef TIMINGS
    timings[3] = LAGraph_WallClockTime ( );
    #endif
    #pragma omp parallel for schedule(static)
    for(i = 0; i < max_deg; ++i)
    {
        rcc[i] = ((double)epd[i]) / ((double)vpd[i] * ((double) vpd[i] - 1.0)) ;
    }
    #ifdef TIMINGS
    timings[4] = LAGraph_WallClockTime ( );
    timings[4] -= timings[3];
    timings[3] -= timings[2];
    timings[2] -= timings[1];
    timings[1] -= timings[0];
    timings[0] -= tic;
    
    print_timings(timings);
    LG_SET_BURBLE(false);
    #endif
    epd = vpd = NULL;
    GRB_TRY (GrB_Vector_new(rccs, GrB_FP64, max_deg));
    GRB_TRY (GxB_Vector_load(
        *rccs, (void **) &rcc, GrB_FP64, max_deg, max_deg * sizeof(double), 
        GrB_DEFAULT, NULL)) ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;    
    #else
    printf("LAGraph_RichClubCoefficient_NoGB needs GB v10\n") ;
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}
#undef TIMINGS
