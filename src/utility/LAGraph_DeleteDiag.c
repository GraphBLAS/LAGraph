//------------------------------------------------------------------------------
// LAGraph_DeleteDiag: removes the diagonal entries from G->A
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

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

    LG_CHECK_INIT (G, msg) ;
    if (G->ndiag == 0)
    {
        // nothing to do
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // delete all properties not affected by the removal of the diagonal
    //--------------------------------------------------------------------------

    LAGraph_BooleanProperty A_structure_is_symmetric = G->A_structure_is_symmetric ;
    LAGraph_TRY (LAGraph_DeleteProperties (G, msg)) ;
    G->A_structure_is_symmetric = A_structure_is_symmetric ;

    //--------------------------------------------------------------------------
    // remove diagonal entries
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_select (G->A, NULL, NULL, GrB_OFFDIAG, G->A, 0, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace, G->ndiag now known to be zero
    //--------------------------------------------------------------------------

    G->ndiag = 0 ;
    return (0) ;
}

