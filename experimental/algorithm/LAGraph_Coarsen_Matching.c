//------------------------------------------------------------------------------
// LAGraph_Coarsen_Matching: Coarsen an undirected graph using an edge matching
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_Coarsen_Matching
(
    // outputs:
    GrB_Matrix *coarsened, // coarsened adjacency
    GrB_Vector *mapping,   // parent mapping (no compression)
    // inputs:
    LAGraph_Graph G,
    int matching_type,     // how to perform the coarsening
    int preserve_mapping,  // preserve original namespace of nodes
    int nlevels,           // #of coarsening levels
    uint64_t seed,         // used for matching
    char *msg
)
{
    LG_CLEAR_MSG ;
    GrB_Matrix A ;                      // adjacency matrix
    GrB_Matrix E ;                      // incidence matrix
    GrB_Matrix E_t ;                    // transpose of incidence
    GrB_Matrix S ;                      // S matrix (S[i][j] -> node j maps to node i in coarsened graph)
    GrB_Matrix S_t ;                    // transpose of S
    GrB_Vector matched_edges ;          // result of maximal matching
    GrB_Vector edge_parent ;            // points to parent (representative) node for each edge
    GrB_Vector node_parent ;            // points to parent (representative) node for each node

    // get A from G
    // run A_to_E to get incidence
    // run matching on incidence
    // result of matching needs to be converted to parent vector
    //      need to take matched edges and arbitrarily choose one endpoint (can use (ANY/MIN)_first_index)
    //      now, we have a mapping from edges to nodes. Need to go from nodes to nodes
    //      multiply with E, use ANY_SECOND to choose the node to point to (at most one entry can match)
    //      now, need to build S matrix from parent vec (write new util function, similar to A_to_E)
    //      now, multiply SAS' to get coarsened adjacency
    //          previously discussed expanding A to EE', but no advantage b/c we already have A

    // check properties (no self-loops, undirected)
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->is_symmetric_structure == LAGraph_TRUE))
    {
        // the structure of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_ASSERT_MSG (false, -105, "G->A must be symmetric") ;
    }
    LG_ASSERT_MSG (G->nself_edges == 0, -107, "G->nself_edges must be zero") ;

    while (nlevels > 0) {
        // get E
        LG_TRY (LAGraph_A_to_E (&E, G, msg)) ;
        GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;
        
        // run maximal matching
        LG_TRY (LAGraph_MaximalMatching (&matched_edges, E, matching_type, seed, msg)) ;

        // make node_parent
        // GRB_TRY (GrB_mxv (node_parent, )

    }

    return (GrB_SUCCESS) ;
}