//------------------------------------------------------------------------------
// LAGraph_Property_RowDegree: determine G->rowdegree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK       \
{                               \
    GrB_free (&S) ;             \
    GrB_free (&x) ;             \
}

#define LAGraph_FREE_ALL        \
{                               \
    LAGraph_FREE_WORK ;         \
    GrB_free (&rowdegree) ;     \
}

#include "LG_internal.h"

int LAGraph_Property_RowDegree  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine G->rowdegree
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix S = NULL ;
    GrB_Vector rowdegree = NULL, x = NULL ;
    LG_CHECK_INIT (G, msg) ;

    if (G->rowdegree != NULL)
    {
        // G->rowdegree already computed
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // determine the size of A
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

    //--------------------------------------------------------------------------
    // compute the rowdegree
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&rowdegree, GrB_INT64, nrows)) ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )

        // x = zeros (ncols,1)
        GrB_TRY (GrB_Vector_new (&x, GrB_INT64, ncols)) ;
        GrB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, ncols, NULL)) ;
        // rowdegree = A*x using the PLUS_PAIR semiring
        GrB_TRY (GrB_mxv (rowdegree, NULL, NULL, GxB_PLUS_PAIR_INT64,
            A, x, NULL)) ;

    #else

        // S<A,struct> = 1
        GrB_TRY (GrB_Matrix_new (&S, GrB_INT64, nrows, ncols)) ;
        GrB_TRY (GrB_Matrix_assign_INT64 (S, A, NULL, (int64_t) 1,
            GrB_ALL, nrows, GrB_ALL, ncols, GrB_DESC_S)) ;
        // rowdegree = reduce (S) to vector, using the PLUS_MONOID
        GrB_TRY (GrB_Matrix_reduce_Monoid (rowdegree, NULL, NULL,
            GrB_PLUS_MONOID_INT64, S, NULL)) ;

    #endif

    G->rowdegree = rowdegree ;
    G->rowdegree_type = GrB_INT64;

    LAGraph_FREE_WORK ;
    return (0) ;
}
