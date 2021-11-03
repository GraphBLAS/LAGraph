//------------------------------------------------------------------------------
// LAGraph_DeleteDiag: removes the diagonal entries from G->A
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK   \
{                           \
    GrB_free (&thunk) ;     \
    GrB_free (&M) ;         \
}

#include "LG_internal.h"

int LAGraph_DeleteDiag      // returns 0 if successful, < 0 if failure
(
    LAGraph_Graph G,        // diagonal entries removed, most properties cleared
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix M = NULL ;
    GrB_Scalar thunk = NULL ;
    LG_CHECK_INIT (G, msg) ;
    if (G->ndiag == 0)
    {
        // nothing to do
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // delete all properties not affected by the removal of the diagonal
    //--------------------------------------------------------------------------

    LAGraph_BooleanProperty A_pattern_is_symmetric = G->A_pattern_is_symmetric ;
    LAGraph_TRY (LAGraph_DeleteProperties (G, msg)) ;
    G->A_pattern_is_symmetric = A_pattern_is_symmetric ;

    //--------------------------------------------------------------------------
    // remove diagonal entries
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;

    GrB_TRY (GrB_Scalar_new (&thunk, GrB_INT64)) ;
    GrB_TRY (GrB_Scalar_setElement (thunk, 0)) ;
    GrB_TRY (GrB_select (A, NULL, NULL, GrB_OFFDIAG, A, thunk, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace, G->ndiag now known to be zero
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    G->ndiag = 0 ;
    return (0) ;
}

