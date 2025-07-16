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
    "void LG_RCC_iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y) \n"\
    "{                                                                      \n"\
        "(*z) = (int64_t)((*x < *y) + (*x <= *y)) ;                         \n"\
    "}"
void LG_RCC_iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y)
{
    (*z) = (int64_t)((*x < *y) + (*x <= *y)) ;
}

#define RICH_CLUB_FORMULA                                                      \
    "void LG_RCC_rich_club_formula(double *z, const int64_t *x, const int64_t *y) \n"\
    "{                                                                      \n"\
    "   (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0)) ; \n"\
    "}"
void LG_RCC_rich_club_formula(double *z, const int64_t *x, const int64_t *y)
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
    int64_t max_deg;
    bool iso = false;

    void *a_space = NULL;
    
    int64_t *node_edges_arr = NULL, *deg_arr = NULL, 
        *epd_arr = NULL, *ones = NULL, 
        *vpd_arr = NULL;
    GrB_Type epd_type = NULL, vpd_type = NULL;
    uint64_t epd_n = 0, vpd_n = 0, epd_size = 0, vpd_size = 0;
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
        &iseq_2lt, (LAGraph_binary_function) (&LG_RCC_iseq_2islt), 
        GrB_INT64, GrB_INT64, GrB_INT64, "LG_RCC_iseq_2islt", ISEQ_2ISLT)) ;
    GRB_TRY (GxB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&LG_RCC_rich_club_formula), 
        GrB_FP64, GrB_INT64, GrB_INT64, 
        "LG_RCC_rich_club_formula", RICH_CLUB_FORMULA)) ;
    #else
    GRB_TRY (GrB_BinaryOp_new(
        &iseq_2lt, (LAGraph_binary_function) (&LG_RCC_iseq_2islt), 
        GrB_INT64, GrB_INT64, GrB_INT64)) ;
    GRB_TRY (GrB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&LG_RCC_rich_club_formula), 
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
    #if LG_SUITESPARSE_GRAPHBLAS_V10
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

        
        GRB_TRY (GrB_Vector_assign_INT64(
            edges_per_deg, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            verts_per_deg, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_Vector_assign_INT64(
            ones_v, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
            
        #ifndef COVERAGE
        GRB_TRY (GrB_Vector_new(&ramp_v, GrB_INT64, edge_vec_nvals + 1)) ;  
        GRB_TRY (GrB_Vector_assign_INT64(
            ramp_v, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_apply (
            ramp_v, NULL, NULL, GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
        #endif

        LG_TRY (LAGraph_FastAssign_Semiring (
            edges_per_deg, NULL, GrB_PLUS_INT64, deg_x, node_edges_x, ramp_v,
            GxB_PLUS_SECOND_INT64, NULL, msg
        )) ;
        LG_TRY (LAGraph_FastAssign_Semiring (
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
        for(int64_t i = max_deg - 1; i > 0; --i)
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
        for(int64_t i = max_deg - 1; i > 0; --i)
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
