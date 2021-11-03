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
    // x = zeros (nrows,1)
    GrB_TRY (GrB_Vector_new (&x, GrB_INT64, nrows)) ;
    GrB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, nrows, NULL)) ;

    if (AT != NULL)
    {
        // G->coldegree = row degree of AT; this will be faster assuming
        // AT is held in a row-oriented format. 
        GrB_TRY (GrB_mxv (coldegree, NULL, NULL, LAGraph_plus_one_int64,
            AT, x, NULL)) ;
    }
    else
    {
        // G->coldegree = column degree of A
        GrB_TRY (GrB_mxv (coldegree, NULL, NULL, LAGraph_plus_one_int64,
            A, x, GrB_DESC_T0)) ;
    }

    G->coldegree = coldegree ;
    G->coldegree_type = GrB_INT64;

    LAGraph_FREE_WORK ;
    return (0) ;
}
