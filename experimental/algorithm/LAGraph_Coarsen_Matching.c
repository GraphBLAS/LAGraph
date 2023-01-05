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
    int combine_weights,   // whether to sum edge weights or just keep the pattern
    int nlevels,           // #of coarsening levels
    uint64_t seed,         // used for matching
    char *msg
)
{
    LG_CLEAR_MSG ;
    LAGraph_Graph G_cpy = NULL ;            // used for the A_to_E function
    GrB_Matrix A = NULL ;                   // adjacency matrix
    GrB_Matrix E = NULL ;                   // incidence matrix
    GrB_Matrix E_t = NULL ;                 // transpose of incidence
    GrB_Matrix S = NULL;                    // S matrix (S[i][j] -> node j maps to node i in coarsened graph)
    GrB_Matrix S_t = NULL ;                 // transpose of S
    GrB_Matrix result = NULL ;              // resulting adjacency matrix
    GrB_Vector matched_edges = NULL ;       // result of maximal matching
    GrB_Vector edge_parent = NULL ;         // points to parent (representative) node for each edge
    GrB_Vector node_parent = NULL ;         // points to parent (representative) node for each node
    GrB_Vector ones = NULL ;                // vector of all 1's

    // get A from G
    // run A_to_E to get incidence
    // run matching on incidence
    // result of matching needs to be converted to parent vector
    //      need to take matched edges and arbitrarily choose one endpoint (can use (ANY/MIN)_first_index)
    //          need to multiply E_t by all ones, mask by chosen edges
    //      now, we have a mapping from edges to nodes. Need to go from nodes to nodes
    //      multiply with E, use ANY_SECOND to choose the node to point to (at most one entry can match)
    //      now, need to build S matrix from parent vec (write new util function, similar to A_to_E)
    //      now, multiply SAS' to get coarsened adjacency
    //          previously discussed expanding A to EE', but no advantage b/c we already have A

    // check properties (no self-loops, undirected)
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        A = G->A ;
    }
    else
    {
        // G is not undirected
        LG_ASSERT_MSG (false, -105, "G must be undirected") ;
    }
    LG_ASSERT_MSG (G->nself_edges == 0, -107, "G->nself_edges must be zero") ;
    
    // make new LAGraph_Graph to use for A_to_E
    LG_TRY (LAGraph_New (&G_cpy, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    LG_TRY (LAGraph_Cached_NSelfEdges (&G_cpy, msg)) ;

    char *A_typename ;
    GrB_Type A_type ;

    // extract type of A to build coarsened matrix with
    LG_TRY (LAGraph_Malloc ((void**) &A_typename, LAGRAPH_MAX_NAME_LEN, sizeof(char), msg)) ;
    LG_TRY (LAGraph_Matrix_TypeName (A_typename, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&A_type, A_typename, msg)) ;
    LG_TRY (LAGraph_Free ((void**) &A_typename, msg)) ;

    GrB_Index num_nodes ;
    GrB_Index num_edges ;

    GRB_TRY (GrB_Matrix_nvals (&num_edges, A)) ;
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;
    num_edges /= 2 ;

    // matching computations are done by casting E entries into FP64 (catch all type)
    GRB_TRY (GrB_Matrix_new (&result, A_type, num_nodes, num_nodes)) ;
    
    GRB_TRY (GrB_Vector_new (&edge_parent, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&node_parent, GrB_UINT64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&ones, GrB_UINT64, num_nodes)) ;

    GRB_TRY (GrB_assign (ones, NULL, NULL, 1, GrB_ALL, num_nodes, NULL)) ;

    while (nlevels > 0) {
        // get E
        LG_TRY (LAGraph_A_to_E (&E, G_cpy, msg)) ;
        GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;
        
        // run maximal matching
        LG_TRY (LAGraph_MaximalMatching (&matched_edges, E, matching_type, seed, msg)) ;

        // make edge_parent
        // want to do E_t * ones and get the first entry for each edge (mask output with matched_edges)
        // now, we have edge_parent (each edge points to its parent node)
        // can do E * edge_parent with min_second to get node_parent

        // make node parent
        // GRB_TRY (GrB_mxv (node_parent, NULL, NULL, min_secondi, A, edge_parent, NULL))
        GRB_TRY (GrB_mxv (node_parent, NULL, NULL, GrB_MIN_SECOND_SEMIRING_UINT64, E, edge_parent, NULL)) ;
        //...
        G_cpy->A = result;
        LG_TRY (LAGraph_Cached_NSelfEdges (&G_cpy, msg)) ;
        nlevels--;
    }
    (*coarsened) = result ;
    (*mapping) = node_parent ;

    LG_TRY (LAGraph_Delete (&G_cpy, msg)) ;

    return (GrB_SUCCESS) ;
}