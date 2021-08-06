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
    #if LG_SUITESPARSE
    GxB_Scalar thunk = NULL ;
    #else
    GrB_Vector thunk = NULL ;       // unused
    #endif
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

    #if LG_SUITESPARSE
    GrB_TRY (GxB_Scalar_new (&thunk, GrB_INT64)) ;
    GrB_TRY (GxB_Scalar_setElement (thunk, 0)) ;
    GrB_TRY (GxB_select (A, NULL, NULL, GxB_OFFDIAG, A, thunk, NULL)) ;
    #else
    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_new (&M, GrB_BOOL, n, n)) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        GrB_TRY (GrB_Matrix_setElement_BOOL (M, 1, i, i)) ;
    }
    // A<!M,struct,replace> = A
    GrB_TRY (GrB_assign (A, M, NULL, A, GrB_ALL, n, GrB_ALL, n, GrB_DESC_RSC)) ;
    #endif

    //--------------------------------------------------------------------------
    // free workspace, G->ndiag now known to be zero
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    G->ndiag = 0 ;
    return (0) ;
}

