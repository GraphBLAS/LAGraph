//------------------------------------------------------------------------------
// LAGraph_Delete: deletes a graph and all its contents
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Delete
(
    // input/output:
    LAGraph_Graph *G,   // the graph to delete; G set to NULL on output.
                        // All internal GrB_Matrix and GrB_Vector objects are
                        // freed, including G->A.  To keep G->A while deleting
                        // the graph G, use:
                        // { A = G->A ; G->A = NULL ; LAGraph_Delete (&G, msg);}
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
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // free the cached contents of the graph
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_DeleteProperties (*G, msg)) ;

    //--------------------------------------------------------------------------
    // delete the primary contents of the graph, and the graph itself
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_free (&((*G)->A))) ;
    LAGraph_Free ((void **) G, NULL) ;
    return (GrB_SUCCESS) ;
}
