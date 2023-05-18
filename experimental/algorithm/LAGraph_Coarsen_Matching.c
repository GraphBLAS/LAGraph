//------------------------------------------------------------------------------
// LAGraph_Coarsen_Matching: Coarsen an undirected graph using an edge matching
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

/*
This method is used to coarsen an undirected graph. The coarsening is based on a maximal matching,
which is handled by LAGraph_MaximalMatching.
The inputs to this algorithm are as follows in order:
1. an LAGraph_Graph containing the target graph to coarsen
2. the type of matching to perform (random, heavy, or light)
3. whether to retain the size of the graph when coarsening. If 1, then nodes that are eliminated by a coarsening step
are turned into singletons. If 0, the nodes are explicitly relabled and the size of the graph is changed.
4. whether edges that are combined during a coarsening step should have their edge weights summed (for an unweighted graph, this
counts the number of combined edges)
5. How many coarsening steps to perform
6. Random seed used for maximal matching
7. msg for LAGraph error reporting

More specifically, the coarsening step involves a reduction from a graph G to G', where we use a bijection f from
nodes in G to nodes in G'. We can consider f(u) to be the parent of node u. 
For each edge (u, v) in G, we add an edge (f(u), f(v)) to G' iff f(u) != f(v). In our case,
this bijection is given by the maximal matching, where for every matched edge, one of the endpoints of the edge is the
parent (representative) of both endpoints, and any node not part of a matched edge is its own parent.  
*/

#include "LG_internal.h"
#include "LAGraphX.h"

#define dbg

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                                \
{                                                   \
    GrB_free(&E) ;                                  \
    GrB_free(&S) ;                                  \
    GrB_free(&matched_edges) ;                      \
    GrB_free(&E_t) ;                                \
    GrB_free(&S_t) ;                                \
    GrB_free(&edge_parent) ;                        \
    GrB_free(&node_parent) ;                        \
    GrB_free(&ones) ;                               \
}                                                   \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    LAGraph_Delete(&G_cpy, msg) ;                   \
    LAGraph_Free((void**) mapping, msg) ;           \
}                                                   \

int LAGraph_Coarsen_Matching
(
    // outputs:
    GrB_Matrix *coarsened,                  // coarsened adjacency
    GrB_Vector **mapping_result,            // array of parent mappings for each level; if preserve_mapping is true, is NULL
    // inputs:
    LAGraph_Graph G,
    LAGraph_Matching_kind matching_type,    // how to perform the coarsening
    int preserve_mapping,                   // preserve original namespace of nodes
    int combine_weights,                    // whether to sum edge weights or just keep the pattern
    GrB_Index nlevels,                      // #of coarsening levels
    uint64_t seed,                          // used for matching
    char *msg
)
{
    LG_CLEAR_MSG ;
    LAGraph_Graph G_cpy = NULL ;            // used for the IncidenceMatrix function
    GrB_Matrix A = NULL ;                   // resulting adjacency matrix (used for output)
    GrB_Matrix E = NULL ;                   // incidence matrix
    GrB_Matrix E_t = NULL ;                 // transpose of incidence
    GrB_Matrix S = NULL ;                   // S matrix (S[i][j] -> node j maps to node i in coarsened graph)
    GrB_Matrix S_t = NULL ;                 // transpose of S
    GrB_Vector matched_edges = NULL ;       // result of maximal matching
    GrB_Vector edge_parent = NULL ;         // points to parent (representative) node for each edge
    GrB_Vector node_parent = NULL ;         // points to parent (representative) node for each node
    GrB_Vector ones = NULL ;                // vector of all 1's

    GrB_Vector *mapping = NULL ;            // resulting array of parent mappings (used for output)

    // check properties (no self-loops, undirected)
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        GRB_TRY (GrB_Matrix_dup (&A, G->A)) ;
    }
    else
    {
        // G is not undirected
        LG_ASSERT_MSG (false, -105, "G must be undirected") ;
    }
    LG_ASSERT_MSG (G->nself_edges == 0, -107, "G->nself_edges must be zero") ;
    
    // make new LAGraph_Graph to use for building incidence matrix
    LG_TRY (LAGraph_New (&G_cpy, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    LG_TRY (LAGraph_Cached_NSelfEdges (G_cpy, msg)) ;

    // set A back
    A = G_cpy->A ;

    GrB_Index num_nodes ;
    GrB_Index num_edges ;

    GRB_TRY (GrB_Matrix_nvals (&num_edges, A)) ;
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;
    num_edges /= 2 ; // since undirected

    GRB_TRY (GrB_Matrix_new (&E_t, GrB_FP64, num_edges, num_nodes)) ;

    GRB_TRY (GrB_Matrix_new (&S_t, GrB_BOOL, num_nodes, num_nodes)) ;

    GRB_TRY (GrB_Vector_new (&edge_parent, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&node_parent, GrB_UINT64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&ones, GrB_UINT64, num_nodes)) ;

    GRB_TRY (GrB_assign (ones, NULL, NULL, 1, GrB_ALL, num_nodes, NULL)) ;

    if (!preserve_mapping) {
        LG_TRY (LAGraph_Malloc ((void**)(&mapping), nlevels, sizeof(GrB_Vector), msg)) ;
    }

    GrB_Index curr_level = 0 ;

    while (nlevels > 0) {
        // get E
        LG_TRY (LAGraph_Incidence_Matrix (&E, G_cpy, msg)) ;
        GRB_TRY (GrB_Matrix_nvals (&num_edges, A)) ;
        num_edges /= 2 ;
        if (!preserve_mapping) {
            GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;
            // resize structures
            // need to resize E_t, edge_parent, node_parent, ones
            GRB_TRY (GrB_Vector_resize (node_parent, num_nodes)) ;
            GRB_TRY (GrB_Vector_resize (ones, num_nodes)) ;
        }
        GRB_TRY (GrB_Matrix_resize (E_t, num_edges, num_nodes)) ;
        GRB_TRY (GrB_Vector_resize (edge_parent, num_edges)) ;

        GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;

        // run maximal matching
        LG_TRY (LAGraph_MaximalMatching (&matched_edges, E, matching_type, seed, msg)) ;

        // make edge_parent
        // want to do E_t * ones and get the first entry for each edge (mask output with matched_edges)
        GRB_TRY (GrB_mxv (edge_parent, matched_edges, NULL, GxB_MIN_SECONDI_INT64, E_t, ones, GrB_DESC_RS)) ;
        // now, we have edge_parent (each edge points to its parent node)
        // can do E * edge_parent with min_second to get node_parent
        GRB_TRY (GrB_mxv (node_parent, NULL, NULL, GrB_MIN_SECOND_SEMIRING_UINT64, E, edge_parent, NULL)) ;
        // need to do eWiseAdd with [1...n] since some nodes might not touch any matched edges, and we want those nodes to point to self
        GRB_TRY (GrB_eWiseAdd (node_parent, node_parent, NULL, GxB_FIRSTI_INT64, ones, node_parent, GrB_DESC_SC)) ;

        #ifdef dbg
            LG_TRY (LAGraph_Vector_Print (node_parent, LAGraph_COMPLETE, stdout, msg)) ;
        #endif
        
        LG_TRY (LAGraph_Parent_to_S (&S, node_parent, preserve_mapping, msg)) ;
        if (!preserve_mapping) {
            // need to resize S_t
            GrB_Index S_rows, S_cols ;
            GRB_TRY (GrB_Matrix_nrows (&S_rows, S)) ;
            GRB_TRY (GrB_Matrix_ncols (&S_cols, S)) ;
            GRB_TRY (GrB_Matrix_resize (S_t, S_cols, S_rows)) ;
            // resize result
            GRB_TRY (GrB_Matrix_resize (A, S_rows, S_rows)) ;
        }
        GRB_TRY (GrB_transpose (S_t, NULL, NULL, S, NULL)) ;

        GrB_Semiring semiring = combine_weights ? GrB_PLUS_TIMES_SEMIRING_FP64 : LAGraph_any_one_bool ;
        GRB_TRY (GrB_mxm (S, NULL, NULL, semiring, S, A, NULL)) ;
        GRB_TRY (GrB_mxm (A, NULL, NULL, semiring, S, S_t, NULL)) ;
        G_cpy->A = A ;
        LG_TRY (LAGraph_Cached_NSelfEdges (G_cpy, msg)) ;
        // parent nodes for matched edges will form self-edges; need to delete
        LG_TRY (LAGraph_DeleteSelfEdges (G_cpy, msg)) ;
        A = G_cpy->A ;
        // want to free before we reassign what they point to
        GRB_TRY (GrB_free (&E)) ;
        GRB_TRY (GrB_free (&S)) ;
        GRB_TRY (GrB_free (&matched_edges)) ;

        // record a deep copy of the current node_parent for the current coarsening level
        GRB_TRY (GrB_Vector_dup ((mapping) + (curr_level * sizeof(GrB_Vector)), node_parent)) ;

        nlevels-- ;
        curr_level++ ;
    }
    (*coarsened) = A ;
    (*mapping_result) = mapping ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}