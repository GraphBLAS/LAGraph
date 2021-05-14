//------------------------------------------------------------------------------
// LAGraph_Property_ColDegree: determine G->coldegree
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
    GrB_free (&coldegree) ;     \
}

#include "LG_internal.h"

int LAGraph_Property_ColDegree  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine G->coldegree
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix S = NULL ;
    GrB_Vector coldegree = NULL, x = NULL ;
    LG_CHECK_INIT (G, msg) ;

    if (G->coldegree != NULL || G->kind == LAGRAPH_ADJACENCY_UNDIRECTED)
    {
        // G->coldegree already computed, or not needed
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // determine the size of A
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Matrix AT = G->AT ;
    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

    //--------------------------------------------------------------------------
    // compute the coldegree
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&coldegree, GrB_INT64, ncols)) ;

    if (AT != NULL)
    {

        //----------------------------------------------------------------------
        // G->coldegree = row degree of AT
        //----------------------------------------------------------------------

        #if LG_SUITESPARSE

            // x = zeros (nrows,1)
            GrB_TRY (GrB_Vector_new (&x, GrB_INT64, nrows)) ;
            GrB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, nrows, NULL)) ;
            // coldegree = AT*x using the PLUS_PAIR semiring
            GrB_TRY (GrB_mxv (coldegree, NULL, NULL, GxB_PLUS_PAIR_INT64,
                AT, x, NULL)) ;

        #else

            // S<AT,struct> = 1
            GrB_TRY (GrB_Matrix_new (&S, GrB_INT64, ncols, nrows)) ;
            GrB_TRY (GrB_Matrix_assign_INT64 (S, AT, NULL, (int64_t) 1,
                GrB_ALL, ncols, GrB_ALL, nrows, GrB_DESC_S)) ;
            // coldegree = reduce (S) to vector, using the PLUS_MONOID
            GrB_TRY (GrB_Matrix_reduce_Monoid (coldegree, NULL, NULL,
                GrB_PLUS_MONOID_INT64, S, NULL)) ;

        #endif

    }
    else
    {

        //----------------------------------------------------------------------
        // G->coldegree = column degree of A
        //----------------------------------------------------------------------

        #if LG_SUITESPARSE

            // x = zeros (nrows,1)
            GrB_TRY (GrB_Vector_new (&x, GrB_INT64, nrows)) ;
            GrB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, nrows, NULL)) ;
            // coldegree = A'*x using the PLUS_PAIR semiring
            GrB_TRY (GrB_mxv (coldegree, NULL, NULL, GxB_PLUS_PAIR_INT64,
                A, x, GrB_DESC_T0)) ;

        #else

            // S<A,struct> = 1
            GrB_TRY (GrB_Matrix_new (&S, GrB_INT64, nrows, ncols)) ;
            GrB_TRY (GrB_Matrix_assign_INT64 (S, A, NULL, (int64_t) 1,
                GrB_ALL, nrows, GrB_ALL, ncols, GrB_DESC_S)) ;
            // coldegree = reduce (S') to vector, using the PLUS_MONOID
            GrB_TRY (GrB_Matrix_reduce_Monoid (coldegree, NULL, NULL,
                GrB_PLUS_MONOID_INT64, S, GrB_DESC_T0)) ;

        #endif

    }

    G->coldegree = coldegree ;
    G->coldegree_type = GrB_INT64;

    LAGraph_FREE_WORK ;
    return (0) ;
}
