//------------------------------------------------------------------------------
// LAGraph_DeleteSelfEdges: removes the diagonal entries from G->A
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DeleteSelfEdges
(
    // input/output:
    LAGraph_Graph G,    // diagonal entries removed, most cached properties cleared
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;
    if (G->nself_edges == 0)
    {
        // nothing to do
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // delete all cached properties not affected by the removal of the diagonal
    //--------------------------------------------------------------------------

    LAGraph_Boolean is_symmetric_structure = G->is_symmetric_structure ;
    LG_TRY (LAGraph_DeleteCached (G, msg)) ;
    G->is_symmetric_structure = is_symmetric_structure ;

    //--------------------------------------------------------------------------
    // remove diagonal entries
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_select (G->A, NULL, NULL, GrB_OFFDIAG, G->A, 0, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace, G->nself_edges now known to be zero
    //--------------------------------------------------------------------------

    G->nself_edges = 0 ;
    return (GrB_SUCCESS) ;
}

