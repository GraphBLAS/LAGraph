//------------------------------------------------------------------------------
// LAGraph_FW: Floyd-Warshall method: all pairs shortest paths
//------------------------------------------------------------------------------

// The input is a square unsymmetric matrix G, for a directed graph.  G can be
// of any type.  If it is real (float or double), a 64-bit integer, or an
// unsigned 32-bit integer, then the output is of type GrB_FP64.  Otherwise,
// the output is of type GrB_INT32.

// TODO consider giving the user control over the data type of D to use.

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
    GrB_Type   *D_type      // type of D
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

    // determine the type of the input and output graphs
    GrB_Type gtype, otype ;
    GrB_BinaryOp op ;
    GrB_UnaryOp idop ;
    GrB_Semiring semiring ;
    GRB_TRY (GxB_Matrix_type (&gtype, G)) ;

    if (gtype == GrB_FP64 || gtype == GrB_FP32 ||
        gtype == GrB_INT64 || gtype == GrB_UINT64 ||
        gtype == GrB_UINT32)
    {
        otype = GrB_FP64 ;
        semiring = GxB_MIN_PLUS_FP64 ;
        op = GrB_MIN_FP64 ;
        idop = GrB_IDENTITY_FP64 ;
    }
    else
    {
        otype = GrB_INT32 ;
        semiring = GxB_MIN_PLUS_INT32 ;
        op = GrB_MIN_INT32 ;
        idop = GrB_IDENTITY_INT32 ;
    }

    GrB_Index n, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, G)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, G)) ;
    LG_ASSERT (n == ncols, GrB_INVALID_VALUE) ;

    (*D_type) = otype ;
    GRB_TRY (GrB_Matrix_new (&D_matrix, otype, n, n)) ;
    GRB_TRY (LG_SET_FORMAT_HINT (D_matrix, LG_BITMAP)) ;
    GRB_TRY (GrB_Matrix_new (&A, otype, n, 1)) ;
    GRB_TRY (GrB_Matrix_new (&B, otype, 1, n)) ;

    GRB_TRY (GrB_apply (D_matrix, GrB_NULL, GrB_NULL, idop, G, GrB_NULL)) ;

    for (GrB_Index k = 0 ; k < n ; k++)
    {
        // A = D(:,k), the kth column
        GRB_TRY (GrB_extract (A, GrB_NULL, GrB_NULL, D_matrix,
            GrB_ALL, n, &k, 1, GrB_NULL)) ;
        // B = D(k,:), the kth row
        GRB_TRY (GrB_extract (B, GrB_NULL, GrB_NULL, D_matrix,
            &k, 1, GrB_ALL, n, GrB_NULL)) ;
        // D = min (D,A*B) with "*" being the min-plus semiring
        GRB_TRY (GrB_mxm (D_matrix, GrB_NULL, op, semiring, A, B, GrB_NULL)) ;
    }

    LG_FREE_WORK ;
    (*D) = D_matrix ;
    return (GrB_SUCCESS) ;
}