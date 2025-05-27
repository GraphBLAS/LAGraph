//------------------------------------------------------------------------------
// LAGraph_MaximalMatching: maximal matching using an adaptation of Luby's MIS algorithm
// on a line graph. Derived from LAGraph_MaximalIndependentSet
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

// TODO: ready for src? need a vanilla non-GxB, and incidence graphs.

/*
Uses a modified version of Luby's MIS algorithm
Major algorithm steps:
- Compute score for each edge
- Find max score neighbor of each edge (*)
- Retain edges with score == max score neighbor (*)
- Add retained edges to result
- Remove retained edges and their neighbors from the graph (*)

(*): these steps involve what can be thought as a "2-hop" process that involves two
GrB_mxv's: the first to go from edges to vertices, and the second from vertices back to edges.
Tying both steps together yields a single BFS-like step in the line graph. A important side effect
of this is that the source edge gets included in the result of this 2-hop step, which cannot be avoided
since we do not compute E'E explicitly.

The input to this method is an incidence matrix E, of size n-by-e where the
undirected graph G has n nodes and e edges.  If the kth edge of G is the edge
(i,j), then the column E(:,k) contains two entries:  E(i,k) and E(j,k), which
have the same value.  If the graph G is weighted, then both E(i,k) and E(j,k)
are equal to the weight of the (i,j) edge.  If G is unweighted, then both are
equal to 1 (and the matrix E is thus iso-valued).

The output is vector 'matching' of size e, where matching(k) is present (and
equal to true) if the kth edge appears in the maximal matching.  If (i,j) is
a matched edge, then no other edges of G that are incident on nodes i and j
appear in the matching.

This method requires O(e) space for an undirected graph with e edges
*/

#include "LG_internal.h"
#include "LAGraphX.h"

// #define dbg

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define OPTIMIZE_PUSH_PULL

#define LG_FREE_WORK                        \
{                                           \
    GrB_free(&score) ;                      \
    GrB_free(&candidates) ;                 \
    GrB_free(&Seed) ;                       \
    GrB_free(&node_degree) ;                \
    GrB_free(&degree) ;                     \
    GrB_free(&max_node_neighbor) ;          \
    GrB_free(&max_neighbor) ;               \
    GrB_free(&new_members) ;                \
    GrB_free(&new_neighbors) ;              \
    GrB_free(&new_members_nodes) ;          \
    GrB_free(&new_members_node_degree) ;    \
    GrB_free(&empty) ;                      \
    GrB_free(&weight) ;                     \
}                                           \

#define LG_FREE_ALL                         \
{                                           \
    LG_FREE_WORK ;                          \
    GrB_free(&result) ;                     \
}                                           \

#define MAX_FAILURES 50

int LAGraph_MaximalMatching
(
    // outputs:
    GrB_Vector *matching,                 // pointer to output vector
    // inputs:
    GrB_Matrix E,                         // incidence
    GrB_Matrix E_t,                       // incidence transposed (if incorrect, results are undefined)
    LAGraph_Matching_kind matching_type,  // 0 (random), 1 (heavy weight matching), 2 (light weight matching)
    uint64_t seed,                        // random number seed
    char *msg
)
{
#if LAGRAPH_SUITESPARSE
    LG_CLEAR_MSG ;

    if ((matching == NULL) || (E == NULL) || (E_t == NULL)) {
        return GrB_NULL_POINTER ;
    }

    GrB_Vector score = NULL ;                   // score for each edge. Computed according to matching_type
    GrB_Vector weight = NULL ;                  // weight of each edge
    GrB_Vector candidates = NULL ;              // set of candidate edges
    GrB_Vector Seed = NULL ;                    // random number seed vector
    GrB_Vector node_degree = NULL ;             // intermediate result for computing edge degree; degree of the node. Only computed once
    GrB_Vector degree = NULL ;                  // edge degree; number of incident edges. Only computed once
    GrB_Vector max_node_neighbor = NULL ;       // intermediate result for computing max edge neighbor; max edge touching a node
    GrB_Vector max_neighbor = NULL ;            // max neighbor of an edge (including itself)
    GrB_Vector new_members = NULL ;             // new edges to include in the matching
    GrB_Vector new_neighbors = NULL ;           // union of new members and their neighbor edges
    GrB_Vector new_members_nodes = NULL ;       // nodes touching an edge in new_members
    GrB_Vector new_members_node_degree = NULL ; // node degree considering only new members
    GrB_Vector result = NULL ;                  // resulting matching    
    GrB_Vector empty = NULL ;                   // empty vector

    GrB_Index num_edges ;
    GrB_Index num_nodes ;

    char typename[LAGRAPH_MAX_NAME_LEN] ;
    GrB_Type type ;
    LG_TRY (LAGraph_Matrix_TypeName (typename, E, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;


    GRB_TRY (GrB_Matrix_nrows (&num_nodes, E)) ;
    GRB_TRY (GrB_Matrix_ncols (&num_edges, E)) ;

    GRB_TRY (GrB_Vector_new (&candidates, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&Seed, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&score, GrB_FP64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&weight, type, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&node_degree, GrB_UINT64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&degree, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&max_node_neighbor, GrB_FP64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&max_neighbor, GrB_FP64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&new_members, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&new_neighbors, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&new_members_nodes, GrB_BOOL, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&new_members_node_degree, GrB_UINT64, num_nodes)) ;
    GRB_TRY (GrB_Vector_new (&result, GrB_BOOL, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&empty, GrB_BOOL, num_edges)) ;
    
    GRB_TRY (GrB_assign (Seed, NULL, NULL, 0, GrB_ALL, num_edges, NULL)) ;

    LG_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;

    // initially all edges are considered
    GRB_TRY (GrB_assign (candidates, NULL, NULL, true, GrB_ALL, 
        num_edges, NULL)) ;
     
    GrB_Index ncandidates ;
    GrB_Index nfailures = 0 ; // counts how many iterations have failed due to invalid matchings

    GRB_TRY (GrB_Vector_nvals(&ncandidates, candidates)) ;

    // for each node, counts incident edges
    GRB_TRY (GrB_mxv (node_degree, NULL, NULL, LAGraph_plus_one_uint64, E, candidates, NULL)) ;

    // for each edge, sums incident edges for each node. Each edge has an excess of 2 degree, but it doesn't matter since
    // we care about relative degree
    GRB_TRY (GrB_mxv (degree, NULL, NULL, LAGraph_plus_second_uint64, E_t, node_degree, NULL)) ;

    GRB_TRY (GrB_reduce (weight, NULL, NULL, GrB_MAX_MONOID_FP64, E_t, NULL)) ; // use ANY ?

    double sparsity_thresh = 
    #ifdef OPTIMIZE_PUSH_PULL
        0.04 ;
    #else
        1.0;
    #endif

    #if defined ( COVERAGE )
        int kount = 0 ;
    #endif

    while (ncandidates > 0) {
        // first just generate the scores again
        GRB_TRY (GrB_eWiseMult (score, candidates, NULL, GrB_DIV_FP64, Seed, degree, GrB_DESC_RS)) ;

        // for light matching, can multiply scores by 1 / (edge weight)
        if (matching_type == LAGraph_Matching_heavy) {
            // heavy
            GRB_TRY (GrB_eWiseMult (score, NULL, NULL, GrB_TIMES_FP64, score, weight, NULL)) ;
        } else if (matching_type == LAGraph_Matching_light) {
            // light
            GRB_TRY (GrB_eWiseMult (score, NULL, NULL, GrB_DIV_FP64, score, weight, NULL)) ;
        }

        // the actual edge selection is common regardless of matching type

        // intermediate result. Max score edge touching each node
        // don't need to clear this out first because we populate the result for all nodes
        if (ncandidates > sparsity_thresh * num_edges) {
            GRB_TRY (LG_SET_FORMAT_HINT (score, LG_BITMAP)) ;
            GRB_TRY (GrB_mxv (max_node_neighbor, NULL, NULL, GrB_MAX_SECOND_SEMIRING_FP64, E, score, NULL)) ;
        } else {
            GRB_TRY (LG_SET_FORMAT_HINT (score, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (max_node_neighbor, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP64, score, E_t, NULL)) ;
        }

        GrB_Index node_nvals ;
        GRB_TRY (GrB_Vector_nvals (&node_nvals, max_node_neighbor)) ;

        // Max edge touching each candidate edge, including itself
        if (node_nvals > sparsity_thresh * num_nodes) {
            GRB_TRY (LG_SET_FORMAT_HINT (max_node_neighbor, LG_BITMAP)) ;
            GRB_TRY (GrB_mxv (max_neighbor, candidates, NULL, GrB_MAX_SECOND_SEMIRING_FP64, E_t, max_node_neighbor, GrB_DESC_RS)) ;
        } else {
            GRB_TRY (LG_SET_FORMAT_HINT (max_node_neighbor, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (max_neighbor, candidates, NULL, GrB_MAX_FIRST_SEMIRING_FP64, max_node_neighbor, E, GrB_DESC_RS)) ;
        }
        // Note that we are using the GE operator and not G, since max_neighbor includes the self score
        // correctness: both score and max_neighbor only have entries for candidates, so no non-candidate members are produced
        // GRB_TRY (GrB_assign (new_members, NULL, NULL, empty, GrB_ALL, num_edges, NULL)) ; // just experimenting
        GRB_TRY (GrB_eWiseAdd (new_members, NULL, NULL, GrB_GE_FP64, score, max_neighbor, NULL)) ;

        // makes new_members structural
        GRB_TRY (GrB_select (new_members, NULL, NULL, GrB_VALUEEQ_BOOL, new_members, true, NULL)) ; 
    #ifdef dbg
        printf("new members for ncandidates = %lld:\n", ncandidates);
        LAGRAPH_TRY (LAGraph_Vector_Print (new_members, LAGraph_SHORT, stdout, msg)) ;
    #endif

        GrB_Index new_members_nvals ;
        GRB_TRY (GrB_Vector_nvals (&new_members_nvals, new_members)) ;

        // check if any node has > 1 edge touching it. 
        if (new_members_nvals > sparsity_thresh * num_edges) {
            GRB_TRY (LG_SET_FORMAT_HINT (new_members, LG_BITMAP)) ;
            GRB_TRY (GrB_mxv (new_members_node_degree, NULL, NULL, LAGraph_plus_one_uint64, E, new_members, NULL)) ;
        } else {
            GRB_TRY (LG_SET_FORMAT_HINT (new_members, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (new_members_node_degree, NULL, NULL, LAGraph_plus_one_uint64, new_members, E_t, NULL)) ;
        }

        GrB_Index max_degree ; 
        GRB_TRY (GrB_reduce (&max_degree, NULL, GrB_MAX_MONOID_UINT64, new_members_node_degree, NULL)) ;

    #if defined ( COVERAGE )
        if (num_nodes == 20 && kount++ == 1) max_degree = 2 ;
        if (num_nodes == 30 && kount++ == 0) max_degree = 2 ;
    #endif

        if (max_degree > 1) {
            nfailures++ ;
            LG_ASSERT_MSG (nfailures <= MAX_FAILURES, LAGRAPH_CONVERGENCE_FAILURE,
                "method has stalled") ;
            // regen seed and seed vector
            LG_TRY (LAGraph_Random_Seed (Seed, seed + nfailures, msg)) ;
            continue ;
        }
        // add new members to result and remove from candidates
        // also want to remove all adjacent edges in new_members from candidates
        GRB_TRY (GrB_assign (result, new_members, NULL, true, GrB_ALL, num_edges, GrB_DESC_S)) ; 
        // to include neighbor edges, need to compute new_neighbors
        // to do this, we need to compute the intermediate result new_members_nodes
        if (new_members_nvals > sparsity_thresh * num_edges) {
            GRB_TRY (LG_SET_FORMAT_HINT (new_members, LG_BITMAP)) ;
            GRB_TRY (GrB_mxv (new_members_nodes, NULL, NULL, LAGraph_any_one_bool, E, new_members, NULL)) ;
        } else {
            GRB_TRY (LG_SET_FORMAT_HINT (new_members, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (new_members_nodes, NULL, NULL, LAGraph_any_one_bool, new_members, E_t, NULL)) ;
        }

        GRB_TRY (GrB_Vector_nvals (&node_nvals, new_members_nodes)) ;

        if (node_nvals > sparsity_thresh * num_nodes) {
            GRB_TRY (LG_SET_FORMAT_HINT (new_members_nodes, LG_BITMAP)) ;
            GRB_TRY (GrB_mxv (new_neighbors, NULL, NULL, LAGraph_any_one_bool, E_t, new_members_nodes, NULL)) ;
        } else {
            GRB_TRY (LG_SET_FORMAT_HINT (new_members_nodes, LG_SPARSE)) ;
            GRB_TRY (GrB_vxm (new_neighbors, NULL, NULL, LAGraph_any_one_bool, new_members_nodes, E, NULL)) ;
        }

    #ifdef dbg
        LAGRAPH_TRY (LAGraph_Vector_Print (new_neighbors, LAGraph_SHORT, stdout, msg)) ;
    #endif
        // removes the union of new_members and their neighbors
        GRB_TRY (GrB_assign (candidates, new_neighbors, NULL, empty, GrB_ALL, num_edges, GrB_DESC_S)) ;

    #ifdef dbg
        printf("candidates:\n");
        LAGRAPH_TRY (LAGraph_Vector_Print (candidates, LAGraph_SHORT, stdout, msg)) ;
    #endif
        GrB_Index last_ncandidates = ncandidates ;

        GrB_Vector_nvals(&ncandidates, candidates) ;
        
        // advance seed vector
        LG_TRY (LAGraph_Random_Next (Seed, msg)) ;

    #if defined ( COVERAGE )
        if (num_nodes == 50 && kount++ == 0)
        {
            // hack the Seed vector
            GRB_TRY (GrB_assign (Seed, NULL, NULL, 42, GrB_ALL, num_edges, NULL)) ;
        }
    #endif
    }

    (*matching) = result ;
    
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
