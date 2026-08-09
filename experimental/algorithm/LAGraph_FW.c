//------------------------------------------------------------------------------
// LAGraph_FW: Floyd-Warshall method: all pairs shortest paths
//------------------------------------------------------------------------------

// The input is a square unsymmetric matrix G, for a directed graph.  G can be
// of any type.  The output is always of type GrB_FP64, regardless of the input
// type, since shortest path distances are naturally floating point and callers
// such as LAGr_ClosenessCentrality require FP64 output.
// TODO: revisit in future; consider giving the user control over the data type of D to use.

// G(i,j) is the edge weight for edge (i,j).  D(i,j) on output is the length of
// the shortest path from node i to j, if the entry is present.  If D(i,j) is
// not present then there is no path from i to j.  The shortest path itself
// is not returned.

// Negative weights are OK, unless there is a negative weight cycle.  In
// that case, the output is undefined.

#define LG_FREE_WORK            \
{                               \
    GrB_free (&A) ;             \
    GrB_free (&B) ;             \
}

#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
    GrB_free (&D_matrix) ;      \
}

#include "LG_internal.h"
#include "LAGraphX.h"

GrB_Info LAGraph_FW
(
    const GrB_Matrix G,     // input graph, with edge weights
    GrB_Matrix *D,          // output graph, created on output
    GrB_Type   *D_type      // type of D (always GrB_FP64)
)
{
    GrB_Info info ;
    char *msg = NULL ;
    GrB_Matrix D_matrix = NULL, A = NULL, B = NULL ;

    // make sure outputs and input are valid
    LG_ASSERT (G != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (D != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (D_type != NULL, GrB_NULL_POINTER) ;
    (*D) = NULL ;
    (*D_type) = NULL ;

    // output is always FP64
    GrB_BinaryOp op       = GrB_MIN_FP64 ;
    GrB_UnaryOp  idop     = GrB_IDENTITY_FP64 ;
    GrB_Semiring semiring = GxB_MIN_PLUS_FP64 ;

    GrB_Index n, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, G)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, G)) ;
    LG_ASSERT (n == ncols, GrB_INVALID_VALUE) ;

    (*D_type) = GrB_FP64 ;
    GRB_TRY (GrB_Matrix_new (&D_matrix, GrB_FP64, n, n)) ;
    GRB_TRY (LG_SET_FORMAT_HINT (D_matrix, LG_BITMAP)) ;
    GRB_TRY (GrB_Matrix_new (&A, GrB_FP64, n, 1)) ;
    GRB_TRY (GrB_Matrix_new (&B, GrB_FP64, 1, n)) ;

    // copy G into D, casting to FP64
    GRB_TRY (GrB_apply (D_matrix, GrB_NULL, GrB_NULL, idop, G, GrB_NULL)) ;

    // Set D(i,i) = 0 for all i. This overrides any self-edge weights in G
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        GRB_TRY (GrB_Matrix_setElement_FP64 (D_matrix, 0.0, i, i)) ;
    }

    for (GrB_Index k = 0 ; k < n ; k++)
    {
        // A = D(:,k), the kth column
        GRB_TRY (GrB_extract (A, GrB_NULL, GrB_NULL, D_matrix,
            GrB_ALL, n, &k, 1, GrB_NULL)) ;
        // B = D(k,:), the kth row
        GRB_TRY (GrB_extract (B, GrB_NULL, GrB_NULL, D_matrix,
            &k, 1, GrB_ALL, n, GrB_NULL)) ;
        // D = min (D, A*B) using the min-plus semiring
        GRB_TRY (GrB_mxm (D_matrix, GrB_NULL, op, semiring, A, B, GrB_NULL)) ;
    }

    LG_FREE_WORK ;
    (*D) = D_matrix ;
    return (GrB_SUCCESS) ;
}