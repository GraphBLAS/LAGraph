//------------------------------------------------------------------------------
// LAGraph_DeleteCached: deletes the cached properties of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DeleteCached
(
    // input/output:
    LAGraph_Graph G,    // G stays valid, only cached properties are freed
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    if (G == NULL)
    {
        // success: nothing to do
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // free all cached properties of the graph
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_free (&(G->AT))) ;
    GRB_TRY (GrB_free (&(G->row_degree))) ;
    GRB_TRY (GrB_free (&(G->col_degree))) ;
    GRB_TRY (GrB_free (&(G->emin))) ;
    GRB_TRY (GrB_free (&(G->emax))) ;

    //--------------------------------------------------------------------------
    // clear the cached scalar properties of the graph
    //--------------------------------------------------------------------------

    G->structure_is_symmetric = LAGRAPH_UNKNOWN ;
    G->emin_state = LAGRAPH_UNKNOWN ;
    G->emax_state = LAGRAPH_UNKNOWN ;
//  G->nonzero = LAGRAPH_UNKNOWN ;
    G->ndiag = LAGRAPH_UNKNOWN ;
    return (GrB_SUCCESS) ;
}
