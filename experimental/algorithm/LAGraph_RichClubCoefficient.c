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

// Contributed by Gabriel A. Gomez, Texas A&M University

//------------------------------------------------------------------------------

// TODO: ready for src

// Get the rich club coefficient of a graph.

// Given a Symetric Graph with no self edges, LAGraph_RichClubCoefficient will
// calculate the rich club coefficients of the graph. 

// The values will be output as a sparse GrB_Vector, the rich club coefficient 
// of k will be found at the closest entry at or above k.

// The G->out_degree cached property must be defined for this method.

// References:

// Julian J. McAuley, Luciano da Fontoura Costa, and Tibério S. Caetano, "The 
// rich-club phenomenon across complex network hierarchies", Applied Physics 
// Letters Vol 91 Issue 8, August 2007. https://arxiv.org/abs/physics/0701290

#include "LG_internal.h"
#include "LAGraphX.h"

void LG_RCC_iseq_2islt(int64_t *z, const int64_t *x, const int64_t *y)
{
    (*z) = (int64_t)((*x < *y) + (*x <= *y)) ;
}

LG_JIT_STRING(
void LG_RCC_rich_club_formula(double *z, const int64_t *x, const int64_t *y)
{
    (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0));
}, LG_RCC_RICH_CLUB_FORMULA)

LG_JIT_STRING(
typedef struct {
    int64_t *arr;
} LG_RCC_context ;
, LG_RCC_CONTEXT)

LG_JIT_STRING(
void LG_RCC_compare_row_col (
    int64_t *z,
    const bool *x,
    GrB_Index ix,
    GrB_Index jx,
    const bool *y,
    GrB_Index iy,
    GrB_Index jy,
    const LG_RCC_context *contx
) {
    int64_t row_deg = contx->arr[ix] ;
    int64_t col_deg = contx->arr[iy] ;
    (*z) = (int64_t)((row_deg < col_deg) + (row_deg <= col_deg)) ;
}
, LG_RCC_COMPARE_ROW_COL)

// define indexbiop for old or vanila GraphBLAS
#if !LG_SUITESPARSE_GRAPHBLAS_V10
typedef GrB_BinaryOp GxB_IndexBinaryOp;
#endif

//------------------------------------------------------------------------------
// LG_RCC_compute_v10: SuiteSparse GraphBLAS v10+ implementation
//------------------------------------------------------------------------------

// Uses GxB_IndexBinaryOp (LG_RCC_compare_row_col) which reads the row and
// column node indices and looks up their degrees from a context array passed
// via a GrB_Scalar.  This avoids building the intermediate D and A_deg
// matrices.  The scatter of per-node edge counts into per-degree buckets uses
// LAGraph_FastAssign_Semiring.

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                    \
{                                                       \
    GrB_free(&degrees) ;                                \
    GrB_free(&deg_x) ;                                  \
    GrB_free(&node_edges) ;                             \
    GrB_free(&node_edges_x) ;                           \
    GrB_free(&ones_v) ;                                 \
    GrB_free(&edges_per_deg) ;                          \
    GrB_free(&verts_per_deg) ;                          \
    GrB_free(&plus_2le) ;                               \
    GrB_free(&rcCalculation) ;                          \
    GrB_free(&node_compare) ;                           \
    GrB_free(&ramp_v) ;                                 \
    GrB_free(&ctx_s) ;                                  \
    GrB_free(&ctx_type) ;                               \
    GrB_free(&node_compare_idxop) ;                     \
    LAGraph_Free((void **) &ctx.arr, NULL) ;            \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
{                                           \
    LG_FREE_WORK ;                          \
    GrB_free (rccs) ;                       \
}

#if LG_SUITESPARSE_GRAPHBLAS_V10
static int LG_RCC_compute_v10
(
    GrB_Vector *rccs,
    LAGraph_Graph G,
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declarations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL ;       // G->A, the adjacency matrix

    // n degrees vector (degrees = out_degree - 1)
    GrB_Vector degrees = NULL, deg_x = NULL ;

    // n x 1: number of edges each node is "responsible" for
    // (smallest-degree end * 2 + same-degree edges, to handle double-counting)
    GrB_Vector node_edges = NULL, node_edges_x = NULL ;

    // max_degree x 1: edge count / vertex count bucketed by degree
    GrB_Vector edges_per_deg = NULL ;
    GrB_Vector verts_per_deg = NULL ;

    // edge_vec_nvals x 1: vector of ones
    GrB_Vector ones_v = NULL ;

    // ramp vector for LAGraph_FastAssign_Semiring
    GrB_Vector ramp_v = NULL ;

    // [+].[LG_RCC_compare_row_col]: compares row vs. col degree via index op
    GrB_Semiring plus_2le = NULL ;

    // 2E_K / (N_k (N_k - 1))
    GrB_BinaryOp rcCalculation = NULL ;

    // index-binary-op wrapper carrying the degree context
    GrB_BinaryOp node_compare = NULL ;
    GxB_IndexBinaryOp node_compare_idxop = NULL ;

    // context: holds pointer to the raw degree array
    LG_RCC_context ctx = {0} ;
    GrB_Scalar ctx_s = NULL ;
    GrB_Type ctx_type = NULL ;

    GrB_Index n ;
    GrB_Index edge_vec_nvals ;
    int64_t max_deg ;

    GrB_Type epd_type = NULL, vpd_type = NULL ;
    uint64_t epd_n = 0, vpd_n = 0, epd_size = 0, vpd_size = 0 ;
    int epd_h = 0, vpd_h = 0 ;
    int64_t *epd_arr = NULL, *vpd_arr = NULL ;

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    GRB_TRY (GrB_Vector_new (&degrees, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_new (&node_edges, GrB_INT64, n)) ;

    // degrees = G->out_degree - 1
    // Fill out degree vector, to target col_scale mxm on graphs
    // with singletons, scalar value irrelevant.
    GRB_TRY (GrB_Vector_assign_INT64 (
        degrees, NULL, NULL, (int64_t) -1, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign (
        degrees, NULL, GrB_PLUS_INT64, G->out_degree, GrB_ALL, 0, NULL)) ;

    GRB_TRY (GxB_Type_new (
        &ctx_type, sizeof(LG_RCC_context), "LG_RCC_context", LG_RCC_CONTEXT)) ;
    GRB_TRY (GxB_IndexBinaryOp_new (
        &node_compare_idxop, (GxB_index_binary_function) (LG_RCC_compare_row_col),
        GrB_INT64, GrB_BOOL, GrB_BOOL, ctx_type, "LG_RCC_compare_row_col",
        LG_RCC_COMPARE_ROW_COL)) ;
    GRB_TRY (GxB_BinaryOp_new(
        &rcCalculation, (GxB_binary_function) (LG_RCC_rich_club_formula),
        GrB_FP64, GrB_INT64, GrB_INT64,
        "LG_RCC_rich_club_formula", LG_RCC_RICH_CLUB_FORMULA)) ;

    GrB_Type type ;
    uint64_t deg_n ;
    uint64_t deg_size ;
    int handling ;
    GRB_TRY (GxB_Vector_unload (
        degrees, (void **) &ctx.arr, &type, &deg_n, &deg_size, &handling, NULL
    )) ;
    GRB_TRY (GrB_Scalar_new (&ctx_s, ctx_type)) ;
    GRB_TRY (GrB_Scalar_setElement_UDT (ctx_s, &ctx)) ;
    GRB_TRY (GxB_BinaryOp_new_IndexOp (
        &node_compare, node_compare_idxop, ctx_s)) ;

    GRB_TRY (GrB_Semiring_new (&plus_2le, GrB_PLUS_MONOID_INT64, node_compare)) ;

    //--------------------------------------------------------------------------
    // Calculations
    //--------------------------------------------------------------------------
    GRB_TRY (GrB_Vector_reduce_INT64(
        &max_deg, NULL, GrB_MAX_MONOID_INT64, G->out_degree, NULL)) ;
    GRB_TRY (GrB_Vector_new (&edges_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY (GrB_Vector_new (&verts_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY (GrB_Vector_new (rccs, GrB_FP64, max_deg)) ;

    GRB_TRY (GrB_free (&degrees)) ;
    GRB_TRY (GrB_Vector_new (&degrees, GrB_BOOL, deg_n)) ;
    GRB_TRY (GrB_Vector_assign_INT64 (
        degrees, NULL, NULL, (bool) 0, GrB_ALL, 0, NULL)) ;

    // Sum the number of edges each node is "responsible" for.
    GRB_TRY (GrB_mxv (
        node_edges, NULL, GrB_PLUS_INT64, plus_2le, A, degrees, NULL)) ;

    // The rest of this is indexing the number of edges and number of nodes at 
    // each degree and then doing a cummulative sum to know the amount of edges 
    // and nodes at degree geq k.
    GRB_TRY (GrB_Vector_nvals (&edge_vec_nvals, node_edges)) ;

    GRB_TRY (GxB_Vector_load (
        degrees, (void **) &ctx.arr, type, deg_n, deg_size, handling, NULL
    )) ;

    if(n == edge_vec_nvals)
    {
        deg_x = degrees ;
        degrees = NULL ;
        node_edges_x = node_edges ;
        node_edges = NULL ;
    }
    else
    {
        GRB_TRY (GrB_Vector_assign(
            degrees, G->out_degree, NULL, degrees, GrB_ALL, 0, GrB_DESC_RS
        )) ;

        GRB_TRY (GrB_Vector_new (&deg_x, GrB_BOOL, 0)) ;
        GRB_TRY (GrB_Vector_new (&node_edges_x, GrB_BOOL, 0)) ;
        GRB_TRY (GxB_Vector_extractTuples_Vector(
            NULL, deg_x, degrees, NULL
        )) ;
        GRB_TRY (GxB_Vector_extractTuples_Vector(
            NULL, node_edges_x, node_edges, NULL
        )) ;
    }
    GRB_TRY (GrB_Vector_nvals (&edge_vec_nvals, node_edges_x))
    GRB_TRY (GrB_Vector_new (&ones_v, GrB_INT64, edge_vec_nvals)) ;

    GRB_TRY (GrB_Vector_assign_INT64(
        edges_per_deg, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign_INT64(
        verts_per_deg, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign_INT64(
        ones_v, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;

#ifndef COVERAGE
    GRB_TRY (GrB_Vector_new (&ramp_v, GrB_INT64, edge_vec_nvals + 1)) ;
    GRB_TRY (GrB_Vector_assign_INT64(
        ramp_v, NULL, NULL, (int64_t) 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_apply (
        ramp_v, NULL, NULL, GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
#endif

    LG_TRY (LAGraph_FastAssign_Semiring (
        edges_per_deg, NULL, GrB_PLUS_INT64, deg_x, node_edges_x, ramp_v,
        LAGraph_plus_second_int64, NULL, msg
    )) ;
    LG_TRY (LAGraph_FastAssign_Semiring (
        verts_per_deg, NULL, GrB_PLUS_INT64, deg_x, ones_v, ramp_v,
        LAGraph_plus_one_int64, NULL, msg
    )) ;

    GRB_TRY (GxB_Vector_unload(
        edges_per_deg, (void **) &epd_arr, &epd_type,
        &epd_n, &epd_size, &epd_h, NULL)) ;
    GRB_TRY (GxB_Vector_unload(
        verts_per_deg, (void **) &vpd_arr, &vpd_type,
        &vpd_n, &vpd_size, &vpd_h, NULL)) ;

    LG_ASSERT (max_deg == vpd_n && max_deg == epd_n, GrB_INVALID_VALUE) ;
    // TODO: should be a GraphBLAS cumulative-sum primitive
    // run a cummulative sum (backwards) on vpd_arr
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

    // Computes the RCC of a matrix
    GRB_TRY(GrB_eWiseMult(
        *rccs, NULL, NULL, rcCalculation, edges_per_deg, verts_per_deg, NULL
    )) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
#endif

//------------------------------------------------------------------------------
// LG_RCC_compute_vanilla: vanilla GraphBLAS / pre-v10 implementation
//------------------------------------------------------------------------------

// Cannot use index binary ops, so materializes row degrees into A_deg = D * A
// (via GrB_Matrix_diag + mxm), where each entry A_deg(i,j) = out_degree(i)-1.
// Uses the value-only binary op iseq_2lt for mxv.  The scatter of per-node
// edge counts into per-degree buckets uses manual C-array extraction followed
// by GrB_Vector_build.

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                    \
{                                                       \
    GrB_free(&D) ;                                      \
    GrB_free(&A_deg) ;                                  \
    GrB_free(&degrees) ;                                \
    GrB_free(&node_edges) ;                             \
    GrB_free(&edges_per_deg) ;                          \
    GrB_free(&verts_per_deg) ;                          \
    GrB_free(&plus_2le) ;                               \
    GrB_free(&rcCalculation) ;                          \
    GrB_free(&node_compare) ;                           \
    LAGraph_Free(&a_space, NULL) ;                      \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
{                                           \
    LG_FREE_WORK ;                          \
    GrB_free (rccs) ;                       \
}

static int LG_RCC_compute_vanilla
(
    GrB_Vector *rccs,
    LAGraph_Graph G,
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declarations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL ;       // G->A, the adjacency matrix

    // n x n: adjacency matrix with A_deg(i,j) = out_degree(i) - 1
    GrB_Matrix A_deg = NULL ;

    // n x n diagonal matrix with D(i,i) = out_degree(i) - 1
    GrB_Matrix D = NULL ;

    // n degrees vector (degrees = out_degree - 1)
    GrB_Vector degrees = NULL ;

    // n x 1: number of edges each node is "responsible" for
    // (smallest-degree end * 2 + same-degree edges, to handle double-counting)
    GrB_Vector node_edges = NULL ;

    // max_degree x 1: edge count / vertex count bucketed by degree
    GrB_Vector edges_per_deg = NULL ;
    GrB_Vector verts_per_deg = NULL ;

    // [+].[iseq_2lt]: compares values directly (row degree stored in A_deg)
    GrB_Semiring plus_2le = NULL ;

    // 2E_K / (N_k (N_k - 1))
    GrB_BinaryOp rcCalculation = NULL ;

    // 2 * (x < y) + (x == y)
    GrB_BinaryOp node_compare = NULL ;

    GrB_Index n ;
    GrB_Index edge_vec_nvals ;
    int64_t max_deg ;
    bool iso = false ;

    void *a_space = NULL ;

    int64_t *node_edges_arr = NULL, *deg_arr = NULL, *epd_arr = NULL,
    *ones = NULL, *vpd_arr = NULL ;
    GrB_Index *epd_index = NULL, *vpd_index = NULL ;

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new (&A_deg, GrB_INT64, n, n)) ;

    GRB_TRY (GrB_Vector_new (&degrees, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_new (&node_edges, GrB_INT64, n)) ;

    // degrees = G->out_degree - 1
    // Fill out degree vector, to target col_scale mxm on graphs
    // with singletons, scalar value irrelevant.
    GRB_TRY (GrB_Vector_assign_INT64 (
        degrees, NULL, NULL, (int64_t) -1, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign (
        degrees, NULL, GrB_PLUS_INT64, G->out_degree, GrB_ALL, 0, NULL)) ;

    GRB_TRY (GrB_BinaryOp_new(
        &node_compare, (GxB_binary_function) (LG_RCC_iseq_2islt),
        GrB_INT64, GrB_INT64, GrB_INT64)) ;
    GRB_TRY (GrB_BinaryOp_new(
        &rcCalculation, (GxB_binary_function) (LG_RCC_rich_club_formula),
        GrB_FP64, GrB_INT64, GrB_INT64 )) ;

    GRB_TRY (GrB_Semiring_new (&plus_2le, GrB_PLUS_MONOID_INT64, node_compare)) ;

    //--------------------------------------------------------------------------
    // Calculations
    //--------------------------------------------------------------------------
    GRB_TRY (GrB_Vector_reduce_INT64(
        &max_deg, NULL, GrB_MAX_MONOID_INT64, G->out_degree, NULL)) ;
    GRB_TRY (GrB_Vector_new (&edges_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY (GrB_Vector_new (&verts_per_deg, GrB_INT64, max_deg)) ;
    GRB_TRY (GrB_Vector_new (rccs, GrB_FP64, max_deg)) ;

    // Each edge in the graph gets the value of the degree of its row node
    GRB_TRY (GrB_Matrix_diag (&D, degrees, 0)) ;
    GRB_TRY (GrB_mxm (
        A_deg, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, D, A, NULL)) ;

    // Sum the number of edges each node is "responsible" for.
    GRB_TRY (GrB_mxv (
        node_edges, NULL, GrB_PLUS_INT64, plus_2le, A_deg, degrees, NULL)) ;

    // The rest of this is indexing the number of edges and number of nodes at 
    // each degree and then doing a cummulative sum to know the amount of edges 
    // and nodes at degree geq k.
    GRB_TRY (GrB_Vector_nvals (&edge_vec_nvals, node_edges)) ;

    LG_TRY (LAGraph_Malloc(
        &a_space, edge_vec_nvals * 3 + max_deg * 4, sizeof(int64_t), NULL
    )) ;
    int64_t *T = a_space ;
    deg_arr = T ;           T += edge_vec_nvals ;
    node_edges_arr = T ;    T += edge_vec_nvals ;
    ones = T ;              T += edge_vec_nvals ;
    epd_arr = T ;           T += max_deg ;
    vpd_arr = T ;           T += max_deg ;
    epd_index = (GrB_Index *) T ;         T += max_deg ;
    vpd_index = (GrB_Index *) T ;         T += max_deg ;

    #pragma omp parallel for schedule(static)
    for(uint64_t i = 0; i < edge_vec_nvals; ++i)
    {
        ones[i] = 1ll ;
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
        edges_per_deg, (GrB_Index *) deg_arr, node_edges_arr, edge_vec_nvals,
        GrB_PLUS_INT64)) ;
    GRB_TRY (GrB_Vector_build_INT64 (
        verts_per_deg, (GrB_Index *) deg_arr, ones, edge_vec_nvals, GrB_PLUS_INT64)) ;
    GRB_TRY (GrB_Vector_assign_INT64(
        edges_per_deg, edges_per_deg, NULL, (int64_t) 0,
        GrB_ALL, 0, GrB_DESC_SC)) ;
    GRB_TRY (GrB_Vector_assign_INT64(
        verts_per_deg, verts_per_deg, NULL, (int64_t) 0,
        GrB_ALL, 0, GrB_DESC_SC)) ;

    // Extract into arrays
    GRB_TRY (GrB_Vector_extractTuples_INT64(
        epd_index, epd_arr, (GrB_Index *) &max_deg, edges_per_deg
    )) ;
    GRB_TRY (GrB_Vector_extractTuples_INT64(
        vpd_index, vpd_arr, (GrB_Index *) &max_deg, verts_per_deg
    )) ;
    // TODO: should be a GraphBLAS cumulative-sum primitive
    // run a cummulative sum (backwards) on vpd_arr
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
    epd_index = NULL ;
    vpd_index = NULL ;
    epd_arr = NULL ;
    vpd_arr = NULL ;

    // Computes the RCC of a matrix
    GRB_TRY(GrB_eWiseMult(
        *rccs, NULL, NULL, rcCalculation, edges_per_deg, verts_per_deg, NULL
    )) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_RichClubCoefficient: public dispatcher
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK { }

#undef  LG_FREE_ALL
#define LG_FREE_ALL  { }

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
    LG_CLEAR_MSG ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rccs != NULL, GrB_NULL_POINTER) ;

    LG_ASSERT_MSG(
        G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE,
        "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,
                   "G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE,
                   "G->nself_edges must be zero") ;

#if LG_SUITESPARSE_GRAPHBLAS_V10
    return LG_RCC_compute_v10 (rccs, G, msg) ;
#else
    return LG_RCC_compute_vanilla(rccs, G, msg) ;
#endif
}
