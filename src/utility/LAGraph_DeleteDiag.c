//------------------------------------------------------------------------------
// LAGraph_DeleteDiag: removes the diagonal entries from G->A
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DeleteDiag
(
    // input/output:
    LAGraph_Graph G,    // diagonal entries removed, most properties cleared
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;
    if (G->ndiag == 0)
    {
        // nothing to do
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // delete all properties not affected by the removal of the diagonal
    //--------------------------------------------------------------------------

    LAGraph_BooleanProperty structure_is_symmetric = G->structure_is_symmetric ;
    LG_TRY (LAGraph_DeleteProperties (G, msg)) ;
    G->structure_is_symmetric = structure_is_symmetric ;

    //--------------------------------------------------------------------------
    // remove diagonal entries
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_select (G->A, NULL, NULL, GrB_OFFDIAG, G->A, 0, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace, G->ndiag now known to be zero
    //--------------------------------------------------------------------------

    G->ndiag = 0 ;
    return (GrB_SUCCESS) ;
}

