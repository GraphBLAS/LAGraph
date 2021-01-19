//------------------------------------------------------------------------------
// LAGraph_Property_ColDegree: determine G->coldegree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_ALL GrB_free (&coldegree) ;

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

    GrB_Vector coldegree = NULL ;
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
    GrB_Index n ;
    GrB_TRY (GrB_Matrix_ncols (&n, A)) ;

    //--------------------------------------------------------------------------
    // compute the coldegree
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&coldegree, GrB_INT64, n)) ;
    GrB_TRY (GrB_assign (coldegree, NULL, NULL, 0, GrB_ALL, n, NULL)) ;

    if (G->AT != NULL)
    {
        // G->coldegree = row degree of AT
        GrB_TRY (GrB_mxv (coldegree, NULL, GrB_PLUS_INT64, GxB_PLUS_PAIR_INT64,
            G->AT, coldegree, NULL)) ;
    }
    else
    {
        GrB_TRY (GrB_mxv (coldegree, NULL, GrB_PLUS_INT64, GxB_PLUS_PAIR_INT64,
            A, coldegree, GrB_DESC_T0)) ;
    }

    G->coldegree = coldegree ;

    return (0) ;
}

