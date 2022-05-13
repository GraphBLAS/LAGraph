//------------------------------------------------------------------------------
// LAGraph_Cached_NSelfEdges: count the # of diagonal entries of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Cached_NSelfEdges
(
    // input/output:
    LAGraph_Graph G,    // graph to compute G->nself_edges
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;

    // already computed
    if (G->nself_edges != LAGRAPH_UNKNOWN)
    {
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // compute G->nself_edges
    //--------------------------------------------------------------------------

    return (LG_nself_edges (&G->nself_edges, G->A, msg)) ;
}

