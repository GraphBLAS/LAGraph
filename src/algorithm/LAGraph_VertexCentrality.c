//------------------------------------------------------------------------------
// LAGraph_VertexCentrality:  centrality dispatch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// This is a Basic algorithm (properties are computed as needed)

// LAGraph_VertexCentrality computes many types of centrality metrics:
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

// TODO: only methods 1 and 2 are currently implemented.
// Neither method considers edge weights, but edge-weights could be handled.

// TODO: do we ensure that sum(centrality) == 1 for all methods?

// TODO: Different metrics may require different input parameters (PageRank
// needs tol, itermax, damping; BC needs # sources or a list of sources, etc).
// It might be hard for a Basic algorithm to pick these parameters by itself,
// and a unified Basic algorithm would need to fit will with all future
// Centrality metrics.  Perhaps we don't write this Basic algorithm yet, and
// add it only when more Centrality metrics are added.

#include "LG_alg_internal.h"

int LAGraph_VertexCentrality    // TODO: write this
(
    // output:
    GrB_Vector *centrality,     // centrality(i): centrality metric of node i
    // input/output:
    LAGraph_Graph G,            // input/output graph
    // input:
    LAGraph_Centrality_Kind kind,    // kind of centrality to compute
//  int accuracy,               // TODO?: 0:quick, 1:better, ... max:exact
    char *msg
)
{
    LG_CLEAR_MSG ;
    return (GrB_NOT_IMPLEMENTED) ;        // TODO
}

