//------------------------------------------------------------------------------
// LAGraph_RichClubCoefficient: rich club coefficient
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Get the rich club coefficient of a graph, also allows for edge randomization
// to normalize the coefficients in a graph.

// Given a Symetric Graph with no self edges, LAGraph_RichClubCoefficient will
// calculate the rich club coefficients of the graph. 

// The values will be output as a sparse GrB_Vector, the rich club coefficient 
// of k will be found at the closeste entry at or above k.

// References:

// Julian J. McAuley, Luciano da Fontoura Costa, and Tibério S. Caetano, “The 
// rich-club phenomenon across complex network hierarchies”, Applied Physics 
// Letters Vol 91 Issue 8, August 2007. https://arxiv.org/abs/physics/0701290

// R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon, “Uniform 
// generation of random graphs with arbitrary degree sequences”, 2006. 
// https://arxiv.org/abs/cond-mat/0312028

#define LG_FREE_WORK                                    \
{                                                       \
    /* free any workspace used here */                  \
    GrB_free(&D) ;                                      \
    GrB_free(&edge_degrees) ;                           \
    GrB_free(&degrees) ;                                \
    GrB_free(&node_edges) ;                             \
    GrB_free(&edges_per_deg) ;                          \
    GrB_free(&verts_per_deg) ;                          \
    GrB_free(&iseq_2lt) ;                               \
    GrB_free(&plus_2le) ;                               \
    GrB_free(&rcCalculation) ;                          \
    LAGraph_Free((void **) &index_edge, NULL) ;         \
    LAGraph_Free((void **) &node_edges_arr, NULL);  \
    LAGraph_Free((void **) &deg_arr, NULL);             \
    LAGraph_Free((void **) &edges_per_deg_arr, NULL);   \
    LAGraph_Free((void **) &ones, NULL);                \
    LAGraph_Free((void **) &deg_vertex_count, NULL);    \
}


#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (rich_club_coefficents) ;      \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

typedef void (*LAGraph_binary_function) (void *, const void *, const void *) ;

/* #define TWO_ONE_ADD                                                         \
    "void two_one_add(int64_t *z, const int64_t *x, const int64_t *y)"   \
    "{"                                                                     \
        "(*z) = 2 * (*x) + (*y) ;"                                          \
    "}"
void two_one_add(int64_t *z, const int64_t *x, const int64_t *y)
{ 
    (*z) = 2 * (*x) + (*y);
}
 */
#define ISEQ_2ISLT                                                          \
    "void iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y)"    \
    "{"                                                                     \
        "(*z) = (int64_t)((*x < *y) + (*x <= *y)) ;"                                   \
    "}"
void iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y)
{
    (*z) = (int64_t)((*x < *y) + (*x <= *y)) ;
}

#define RICH_CLUB_FORMULA                                                      \
    "void rich_club_formula(double *z, const int64_t *x, const int64_t *y)"  \
    "{"                                                                        \
        "(*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0)) ;"  \
    "}"
void rich_club_formula(double *z, const int64_t *x, const int64_t *y)
{
    (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0));
} 

int LAGraph_RichClubCoefficient
(
    // output:
    //rich_club_coefficents(i): rich club coefficents of i
    GrB_Vector *rich_club_coefficents,    

    // input: 
    LAGraph_Graph G, //input graph
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    LG_CLEAR_MSG ;
    //Matrix containing every edge 
    //With an entry cooresponding to the degree of its column
    GrB_Matrix edge_degrees = NULL;

    //A matrix with diagonal entries corresponding to degrees.
    GrB_Matrix D = NULL;

    //degrees of nodes.
    GrB_Vector degrees = NULL;

    // contains the number of edges for which the ith node is
    // the smallest degree node * 2 + # edges w/ same degree as the other node
    // to account for double counting of edges w/ same degree as the other node.
    GrB_Vector node_edges = NULL;

    // the ith entry contains the number of edges whose lowest degree is i.
    GrB_Vector edges_per_deg = NULL;

    // the ith entry contains the number of verticies whose degree is i.
    GrB_Vector verts_per_deg = NULL;

    // 2 * (x < y) + (x == y)
    GrB_BinaryOp iseq_2lt = NULL;

    // [+].[iseq_2lt]
    GrB_Semiring plus_2le;

    // 2E_K / (N_k (N_k -1))
    GrB_BinaryOp rcCalculation = NULL;

    GrB_Matrix A ; // G->A, the adjacency matrix
    GrB_Index n ;
    GrB_Index vi_size ;
    GrB_Index vx_size ;
    GrB_Index edge_vec_nvals;
    GrB_Index deg_vec_size;
    bool iso = false;
    
    int64_t *node_edges_arr = NULL, *deg_arr = NULL, 
        *edges_per_deg_arr = NULL, *ones = NULL, 
        *deg_vertex_count = NULL;
    GrB_Index *index_edge = NULL;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rich_club_coefficents != NULL, GrB_NULL_POINTER);

    //TODO: double check this
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
    GRB_TRY(GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(&edge_degrees, GrB_INT64,n,n)) ;


    GRB_TRY(GrB_Vector_new(&degrees, GrB_INT64, n)) ;
    GRB_TRY(GrB_Vector_new(&node_edges, GrB_INT64, n)) ;

    GRB_TRY(GxB_BinaryOp_new(
        &iseq_2lt, (LAGraph_binary_function) (&iseq_2islt), 
        GrB_INT64, GrB_INT64, GrB_INT64, "iseq_2islt", ISEQ_2ISLT)) ;
    GRB_TRY(GrB_Semiring_new(&plus_2le, GrB_PLUS_MONOID_INT64, iseq_2lt)) ;
    GRB_TRY(GxB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&rich_club_formula), 
        GrB_FP64, GrB_INT64, GrB_INT64, 
        "rich_club_formula", RICH_CLUB_FORMULA)) ;

    // degrees = G->out_degree
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64(
        degrees, NULL, NULL, GrB_MINUS_INT64, G->out_degree, 1, NULL)) ;

    
    // Fill out degree vector, to activate col_scale and row_scale on graphs 
    // with singletons, scalar irrelevant
    GRB_TRY (GrB_Vector_assign_INT64(
        degrees, degrees, NULL, (int64_t) -1, GrB_ALL, 0, GrB_DESC_SC)) ;
    GRB_TRY (GrB_Matrix_diag(&D, degrees, 0)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64(
        degrees, NULL, NULL, GrB_MINUS_INT64, G->out_degree, 1, GrB_DESC_R)) ;
    int64_t max_deg;
    GRB_TRY(GrB_Vector_reduce_INT64(
        &max_deg, NULL, GrB_MAX_MONOID_INT64, G->out_degree, NULL)) ;
    GRB_TRY(GrB_Vector_new(&edges_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY(GrB_Vector_new(&verts_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY(GrB_Vector_new(rich_club_coefficents, GrB_FP64, max_deg)) ;
    //--------------------------------------------------------------------------
    // Calculating time
    //--------------------------------------------------------------------------

    // Each edge in the graph gets the value of the degree of its column node
    GRB_TRY (GrB_mxm(
        edge_degrees, NULL, NULL, GxB_ANY_SECOND_INT64, A, D, NULL)) ;

    // If the nodes of an edge have different degrees, the edge is counted once.
    // If they have the same degree, that edge is double counted. So, we adjust:
    GRB_TRY(GrB_mxm(
        edge_degrees, NULL, NULL, plus_2le, D, edge_degrees, NULL)) ;

    // Sum up the number of edges each node is "responsible" for.
    GRB_TRY(GrB_Matrix_reduce_Monoid(
        node_edges, NULL, NULL, GrB_PLUS_MONOID_INT64, edge_degrees, NULL)) ;

    // The rest of this is indexing the number of edges and number of nodes at 
    // each degree and then doing a cummulative sum to know the amount of edges 
    // and nodes at degree geq k.
    GRB_TRY(GrB_Vector_nvals (&edge_vec_nvals, node_edges)) ;

    // Grab the index and edge count arrays from GBLAS
    // Jumbled NULL so must return sorted. Needed because arrays with 
    // # of edges and # of degrees should line up.
    if(n == edge_vec_nvals)
    {
        GRB_TRY (GxB_Vector_unpack_Full (
            node_edges, (void **)&node_edges_arr, &vx_size, &iso, NULL)) ;
        GRB_TRY (GxB_Vector_unpack_Full (
            degrees, (void **)&deg_arr, &vx_size, &iso, NULL)) ;

        deg_vec_size = n;
    }
    else
    {
        GRB_TRY(GxB_Vector_unpack_CSC(
            node_edges, &index_edge, (void **) &node_edges_arr,
            &vi_size,&vx_size,&iso,&edge_vec_nvals, NULL, NULL)) ;
        LG_TRY(LAGraph_Free((void **)&index_edge, msg)) ;

        GRB_TRY(GxB_Vector_unpack_CSC(
            degrees, &index_edge, (void **) &deg_arr,
            &vi_size,&vx_size,&iso,&deg_vec_size, NULL, NULL)) ; 
        LG_TRY(LAGraph_Free((void **)&index_edge, msg)) ;
    }

    // Build with degrees as indecies and handle duplicates via adition
    GRB_TRY(GrB_Vector_build (
        edges_per_deg, deg_arr, node_edges_arr, deg_vec_size, 
        GrB_PLUS_INT64)) ;
    // TODO: Make ones array in a better way
    LG_TRY(
        LAGraph_Malloc((void **) &ones, deg_vec_size, sizeof(int64_t), NULL)) ;
    for(uint64_t i = 0; i < deg_vec_size; ++i)
        ones[i] = 1;

    GRB_TRY(GrB_Vector_build (
        verts_per_deg, deg_arr, ones, deg_vec_size, GrB_PLUS_INT64)) ;

    GrB_Index *epd_index = NULL,  *vpd_index = NULL;
    GRB_TRY(GxB_Vector_unpack_CSC(
        edges_per_deg, &epd_index, (void **)&edges_per_deg_arr,
        &vi_size, &vx_size, &iso, &edge_vec_nvals, NULL, NULL)) ;
    GRB_TRY(GxB_Vector_unpack_CSC(
        verts_per_deg, &vpd_index, (void **)&deg_vertex_count,
        &vi_size, &vx_size, &iso, &deg_vec_size, NULL, NULL)) ;

    //TODO: parralelize these sums:
    //run a cummulative sum (backwards) on deg_vertex_count
    for(uint64_t i = deg_vec_size - 1; i > 0; --i)
    {
        deg_vertex_count[i-1]+=deg_vertex_count[i];
    }

    //run a cummulative sum (backwards) on edges_per_deg_arr
    for(uint64_t i = edge_vec_nvals - 1; i > 0; --i)
    {
        edges_per_deg_arr[i-1]+=edges_per_deg_arr[i];
    }

    GRB_TRY(GxB_Vector_pack_CSC(
        edges_per_deg, &epd_index, (void **)&edges_per_deg_arr,
        vi_size, vx_size, false, edge_vec_nvals, NULL, NULL));
    GRB_TRY(GxB_Vector_pack_CSC(
        verts_per_deg, &vpd_index, (void **)&deg_vertex_count,
        vi_size, vx_size, false, deg_vec_size, NULL, NULL));

    //Computes the RCC of a matrix
    GRB_TRY(GrB_eWiseMult(*rich_club_coefficents, NULL, NULL, rcCalculation, 
        edges_per_deg, verts_per_deg, NULL)) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
