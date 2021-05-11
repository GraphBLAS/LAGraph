//------------------------------------------------------------------------------
// LAGraph_Property_RowDegree: determine G->rowdegree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_ALL GrB_free (&rowdegree) ;

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

    GrB_Vector rowdegree = NULL ;
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
    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;

    //--------------------------------------------------------------------------
    // compute the rowdegree
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&rowdegree, GrB_INT64, n)) ;
    GrB_TRY (GrB_assign (rowdegree, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_mxv (rowdegree, NULL, GrB_PLUS_INT64,
        GxB_PLUS_PAIR_INT64,        // FIXME
        A, rowdegree, NULL)) ;
    G->rowdegree = rowdegree ;
    G->rowdegree_type = GrB_INT64;

    return (0) ;
}
