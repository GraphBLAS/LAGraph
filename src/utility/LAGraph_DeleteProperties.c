//------------------------------------------------------------------------------
// LAGraph_DeleteProperties: deletes the cached properties of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// TODO: this should also probably set A_pattern_is_symmetric to unknown,
// and ndiag to -1

#include "LG_internal.h"

int LAGraph_DeleteProperties    // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // G stays valid, only cached properties are freed
    char *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    if (G == NULL)
    {
        // success: nothing to do   TODO: or is this an error?
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // free all cached properties of the graph
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_free (&(G->AT))) ;
    GrB_TRY (GrB_free (&(G->rowdegree))) ;
    G->rowdegree_type = NULL;
    GrB_TRY (GrB_free (&(G->coldegree))) ;
    G->coldegree_type = NULL;

    // success
    return (0) ;
}
