//------------------------------------------------------------------------------
// LAGraph_BreadthFirstSearch:  breadth-first search dispatch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LG_alg_internal.h"

// TODO the following is a draft:
int LAGraph_VertexCentrality    // TODO: write this
(
    // output:
    GrB_Vector *centrality,     // centrality(i): centrality metric of node i
    // input/output:
    LAGraph_Graph G,            // input graph
    // input:
    LAGraph_Centrality_Kind kind,    // kind of centrality to compute
//  int accuracy,               // TODO?: 0:quick, 1:better, ... max:exact
    char *msg
)
{
    LG_CLEAR_MSG ;
    return (GrB_NOT_IMPLEMENTED) ;        // TODO
}
