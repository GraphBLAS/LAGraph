//------------------------------------------------------------------------------
// LAGraph_FW: Floyd-Warshall method: all pairs shortest paths
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

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

#define FW_FREE_WORK            \
    GrB_free (&A);              \
    GrB_free (&B);

#define LAGraph_FREE_ALL        \
    FW_FREE_WORK ;              \
    GrB_free (D);

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
GrB_Info LAGraph_FW
(
    const GrB_Matrix G,     // input graph, with edge weights
    GrB_Matrix *D,          // output graph, created on output
    GrB_Type   *D_type      // output type
)
{
#if !defined(LG_SUITESPARSE)
    // GxB_type, semirings and select required
    // FIXME: this can be pure GrB
    return GrB_PANIC;
#else

    GrB_Info info;
    GrB_Matrix A = NULL, B = NULL ;

    // make sure D is valid
    if (D == NULL || D_type == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }
    (*D) = NULL ;

    // determine the type of the input and output graphs
    GrB_Type gtype, otype ;
    GrB_BinaryOp op ;
    GrB_UnaryOp idop ;
    GrB_Semiring semiring ;
    LAGRAPH_OK (GxB_Matrix_type (&gtype, G)) ;
    if (gtype == GrB_FP64 || gtype == GrB_FP32 ||
        gtype == GrB_INT64 || gtype == GrB_UINT64 ||
        gtype == GrB_UINT32)
    {
        // TODO test this
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

    GrB_Index n, ncols;
    LAGRAPH_OK (GrB_Matrix_nrows(&n, G));
    LAGRAPH_OK (GrB_Matrix_ncols(&ncols, G));
    if (n != ncols)
    {
        return (GrB_INVALID_VALUE) ;
    }

    LAGRAPH_OK (GrB_Matrix_new (D,  otype, n, n)) ;
    *D_type = otype;
    LAGRAPH_OK (GrB_Matrix_new (&A, otype, n, 1));
    LAGRAPH_OK (GrB_Matrix_new (&B, otype, 1, n));

    // D = G, with possible typecasting
//  LAGRAPH_OK (GrB_assign (*D, GrB_NULL, GrB_NULL, G, GrB_ALL, n, GrB_ALL, n, GrB_NULL)) ;

// GrB_Matrix_apply (C,Mask,acc,op,A,d)  // C<Mask> = accum (C, op(A))
    LAGRAPH_OK (GrB_apply (*D, GrB_NULL, GrB_NULL, idop, G, GrB_NULL)) ;

    for (GrB_Index i = 0; i < n; i++)
    {
        // A = D (:,i), the ith column
        LAGRAPH_OK (GrB_extract(A, GrB_NULL, GrB_NULL, *D, GrB_ALL, n, &i, 1, GrB_NULL));
        // B = D (i,:), the ith row
        LAGRAPH_OK (GrB_extract(B, GrB_NULL, GrB_NULL, *D, &i, 1, GrB_ALL, n, GrB_NULL));
        // D = min (D,A*B) with "*" being the min-plus semiring
        LAGRAPH_OK (GrB_mxm(*D, GrB_NULL, op, semiring, A, B, GrB_NULL));
    }

    FW_FREE_WORK ;
    return GrB_SUCCESS;
#endif
}
