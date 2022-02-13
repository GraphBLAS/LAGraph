//------------------------------------------------------------------------------
// LAGraph_VertexCentrality:  centrality dispatch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// LAGraph_VertexCentrality is a Basic method for computing many types of
// centrality metrics:
//
//  (1) betweenness centrality
//  (2) pagerank
//  (3) degree (for undirected graphs) or indegree / outdegree (directed graphs)
//  (4) closeness (undirected) or incloseness / outcloseness (directed graphs)
//  (5) eigenvector
//  (6) hubs and authorities
//
// A centrality metric measures the importance of a node in a graph, where on
// output, centrality(i) is the measure of the node's importance based on the
// connectivity of the graph and (optionally) the edge weights of the graph.
// For some methods, exact metrics are very costly to compute for large graphs,
// and thus the methods can be used to compute approximate metrics instead,
// which is much faster.

// TODO: only methods 1 and 2 are currently implemented.  Method (2) currently
// uses the GAP algorithm, which ignores dangling nodes (we need a normal
// pagerank).  No method considers edge weights.

// TODO: do we ensure that sum(centrality) == 1 for all methods?

#include "LG_alg_internal.h"

// TODO the following is a draft:
int LAGraph_VertexCentrality    // TODO: write this
(
    // output:
    GrB_Vector *centrality,     // centrality(i): centrality metric of node i
    // input:
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
