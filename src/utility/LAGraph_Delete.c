//------------------------------------------------------------------------------
// LAGraph_Delete: deletes a graph and all its contents
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Delete      // returns 0 if successful, -1 if failure
(
    LAGraph_Graph *G,   // the graph to delete; G set to NULL on output
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    if (G == NULL || (*G) == NULL)
    {
        // success: nothing to do
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // free the cached contents of the graph
    //--------------------------------------------------------------------------

    LAGraph_TRY (LAGraph_DeleteProperties (*G, msg)) ;

    //--------------------------------------------------------------------------
    // delete the primary contents of the graph, and the graph itself
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_free (&((*G)->A))) ;
    LAGraph_Free ((void **) G) ;

    // success
    return (0) ;
}
