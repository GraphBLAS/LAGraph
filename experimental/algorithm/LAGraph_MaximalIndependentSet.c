//------------------------------------------------------------------------------
// LAGraph_MaximalIndependentSet: maximal independent set, with constraints
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University, ...

//------------------------------------------------------------------------------

// Modified from the GraphBLAS C API Specification, by Aydin Buluc, Timothy
// Mattson, Scott McMillan, Jose' Moreira, Carl Yang.  Based on "GraphBLAS
// Mathematics" by Jeremy Kepner.

#define LAGraph_FREE_WORK           \
{                                   \
    GrB_free (&neighbor_max) ;      \
    GrB_free (&new_members) ;       \
    GrB_free (&new_neighbors) ;     \
    GrB_free (&candidates) ;        \
    GrB_free (&empty) ;             \
    GrB_free (&Seed) ;              \
    GrB_free (&score) ;             \
    GrB_free (&score_zero) ;        \
}

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_FREE_WORK ;             \
    GrB_free (&iset) ;              \
}

#include "LG_internal.h"
#include "LAGraphX.h"

// A variant of Luby's randomized algorithm [Luby 1985]. 

// Given a numeric n x n adjacency matrix A of an unweighted and undirected
// graph (where the value true represents an edge), compute a maximal set of
// independent nodes and return it in a boolean n-vector, mis where
// mis[i] == true implies node i is a member of the set.

// The graph cannot have any self edges, and it must be symmetric.  Self-edges
// (diagonal entries) will cause the method to stall, and thus G->ndiag must be
// zero on input.  G->rowdegree must be present on input.  It must not contain
// any explicit zeros (this is handled by LAGraph_Property_RowDegree).

// Singletons require special treatment.  Since they have no neighbors, their
// score is never greater than the max of their neighbors, so they never get
// selected and cause the method to stall.  To avoid this case they are removed
// from the candidate set at the begining, and added to the independent set.

int LAGraph_MaximalIndependentSet       // maximal independent set
(
    // outputs:
    GrB_Vector *mis,            // mis(i) = true if i is in the set
    // inputs:
    LAGraph_Graph G,            // input graph
    int64_t seed,               // random number seed
    GrB_Vector ignore_node,     // if NULL, no nodes are ignored.  Otherwise
                                // ignore_node(i) = true if node i is to be
                                // ignored, and not treated as a candidate
                                // added to maximal independent set.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector iset = NULL ;            // independent set (output vector)
    GrB_Vector score = NULL ;           // random score for each node
    GrB_Vector score_zero = NULL ;      // where score is zero
    GrB_Vector neighbor_max = NULL ;    // value of max neighbor score
    GrB_Vector new_members = NULL ;     // set of new members to add to iset
    GrB_Vector new_neighbors = NULL ;   // new neighbors to new iset members
    GrB_Vector candidates = NULL ;      // candidate nodes
    GrB_Vector empty = NULL ;           // an empty vector
    GrB_Vector Seed = NULL ;            // random number seed vector
    GrB_Vector degree = NULL ;          // G->rowdegree
    GrB_Matrix A ;                      // G->A, the adjacency matrix
    GrB_Index n ;                       // # of nodes

    LG_CHECK (LAGraph_CheckGraph (G, msg), -102, "graph is invalid") ;
    LG_CHECK (mis == NULL, -103, "mis is null") ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_structure_is_symmetric == LAGRAPH_TRUE))
    {
        // the structure of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -105, "G->A must be symmetric") ;
    }

    degree = G->rowdegree ;
    LG_CHECK (degree == NULL, -106, "G->rowdegree must be defined") ;
    LG_CHECK (G->ndiag != 0, -107, "G->ndiag must be zero") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Vector_new (&neighbor_max, GrB_INT64, n)) ;
    GrB_TRY (GrB_Vector_new (&new_members, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&new_neighbors, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&candidates, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&empty, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&Seed, GrB_INT64, n)) ;
    GrB_TRY (GrB_Vector_new (&score, GrB_INT64, n)) ;
    GrB_TRY (GrB_Vector_new (&score_zero, GrB_INT64, n)) ;
    GrB_TRY (GrB_Vector_new (&iset, GrB_BOOL, n)) ;

    //--------------------------------------------------------------------------
    // remove singletons (nodes of degree zero) and handle ignore_node
    //--------------------------------------------------------------------------

    GrB_Index nonsingletons = 0 ;
    GrB_TRY (GrB_Vector_nvals (&nonsingletons, degree)) ;
    if (nonsingletons == n)
    { 
        if (ignore_node == NULL)
        {
            // all nodes have degree 1 or more; all nodes are candidates
            // candidates (0:n-1) = true
            GrB_TRY (GrB_assign (candidates, NULL, NULL, (bool) true, GrB_ALL,
                n, NULL)) ;
            // Seed vector starts out dense
            // Seed (0:n-1) = 0
            GrB_TRY (GrB_assign (Seed, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
        }
        else
        {
            // all nodes have degree 1 or more, but some nodes are to be
            // ignored.  Use ignore_node as a valued mask.
            // candidates<!ignore_node> = true
            GrB_TRY (GrB_assign (candidates, ignore_node, NULL, (bool) true,
                GrB_ALL, n, GrB_DESC_C)) ;
            // Seed vector starts out sparse
            // Seed{candidates} = 0
            GrB_TRY (GrB_assign (Seed, candidates, NULL, (int64_t) 0, GrB_ALL,
                n, GrB_DESC_S)) ;
        }
    }
    else
    {
        // one or more singleton is present.
        // candidates{degree} = 1
        GrB_TRY (GrB_assign (candidates, degree, NULL, (bool) true,
            GrB_ALL, n, GrB_DESC_S)) ;
        // add all singletons to iset
        // iset{!degree} = 1
        GrB_TRY (GrB_assign (iset, degree, NULL, (bool) true, GrB_ALL, n,
            GrB_DESC_SC)) ;
        if (ignore_node != NULL)
        {
            // one or more singletons are present, and some nodes are to be
            // ignored.  The candidates are all those nodes with degree > 0
            // for which ignore_node(i) is false (or not present).  Delete
            // any candidate i for which ignore_node(i) is true.  Use
            // ignore_node as a valued mask.
            // candidates<ignore_node> = empty
            GrB_TRY (GrB_assign (candidates, ignore_node, NULL, empty,
                GrB_ALL, n, NULL)) ;
            // Delete any ignored nodes from iset
            // iset<ignore_node> = empty
            GrB_TRY (GrB_assign (iset, ignore_node, NULL, empty,
                GrB_ALL, n, NULL)) ;
        }
        // Seed vector starts out sparse
        // Seed{candidates} = 0
        GrB_TRY (GrB_assign (Seed, candidates, NULL, (int64_t) 0, GrB_ALL, n,
            GrB_DESC_S)) ;
    }

    // create the random number seeds
    LAGraph_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;

    //--------------------------------------------------------------------------
    // iterate while there are candidates to check
    //--------------------------------------------------------------------------

    int nstall = 0 ;
    GrB_Index ncandidates ;
    GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
    GrB_Index last_ncandidates = ncandidates ;
    GrB_Index n1 = (GrB_Index) (0.04 * (double) n) ;
    GrB_Index n2 = (GrB_Index) (0.10 * (double) n) ;

    while (ncandidates > 0)
    {
        // compute the score for each node; scale the Seed by the degree
        // score = Seed / degree
        GrB_TRY (GrB_eWiseMult (score, NULL, NULL, GrB_DIV_INT64, Seed, degree,
            NULL)) ;

        // find any score equal to zero, and set it back to the Seed value
        GrB_TRY (GrB_select (score_zero, NULL, NULL, GrB_VALUEEQ_INT64, score,
            (int64_t) 0, NULL)) ;

        GrB_Index nzero ;
        GrB_TRY (GrB_Vector_nvals (&nzero, score_zero)) ;
        if (nzero > 0)
        {
            // Reset any score equal to zero back to its unscaled seed value.
            // This case is very rare and hard to test.
            // score{score == 0} = Seed
            GrB_TRY (GrB_assign (score, score_zero, NULL, Seed, GrB_ALL, n,
                GrB_DESC_S)) ;
        }

        // compute the max score of all candidate neighbors (only candidates
        // have a score, so non-candidate neighbors are excluded)
        // neighbor_max{candidates,replace} = score * A
        GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
        if (ncandidates < n1)
        {
            // push
            // neighbor_max'{candidates,replace} = score' * A
            GrB_TRY (GrB_vxm (neighbor_max, candidates, NULL,
                GrB_MAX_FIRST_SEMIRING_INT64, score, A, GrB_DESC_RS)) ;
        }
        else
        {
            // pull
            // neighbor_max{candidates,replace} = A * score
            GrB_TRY (GrB_mxv (neighbor_max, candidates, NULL,
                GrB_MAX_SECOND_SEMIRING_INT64, A, score, GrB_DESC_RS)) ;
        }

        // select node if its score is > than all its active neighbors
        // new_members = (score > neighbor_max) using set union so that nodes
        // with no neighbors fall through to the output, as true (since no
        // score is equal to zero).
        GrB_TRY (GrB_eWiseAdd (new_members, NULL, NULL, GrB_GT_INT64,
            score, neighbor_max, NULL)) ;

        // drop explicit zeros from new_members
        GrB_TRY (GrB_select (new_members, NULL, NULL, GrB_VALUEEQ_BOOL,
            new_members, (bool) true, NULL)) ;

        // add new members to independent set
        // iset{new_members} = true
        GrB_TRY (GrB_assign (iset, new_members, NULL, (bool) true,
            GrB_ALL, n, GrB_DESC_S)) ;

        // remove new members from set of candidates
        // candidates{new_members} = empty
        GrB_TRY (GrB_assign (candidates, new_members, NULL, empty,
            GrB_ALL, n, GrB_DESC_S)) ;

        // early exit if candidates is empty
        GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
        if (ncandidates == 0) { break ; }

        // Neighbors of new members can also be removed from candidates
        // new_neighbors{candidates,replace} = new_members * A
        GrB_Index n_new_members ;
        GrB_TRY (GrB_Vector_nvals (&n_new_members, new_members)) ;
        if (n_new_members < n2)
        {
            // push
            // new_neighbors{candidates,replace} = new_members' * A
            GrB_TRY (GrB_vxm (new_neighbors, candidates, NULL,
                LAGraph_structural_bool, new_members, A, GrB_DESC_RS)) ;
        }
        else
        {
            // pull
            // new_neighbors{candidates,replace} = A * new_members
            GrB_TRY (GrB_mxv (new_neighbors, candidates, NULL,
                LAGraph_structural_bool, A, new_members, GrB_DESC_RS)) ;
        }

        // remove new neighbors of new members from set of candidates
        // candidates{new_neighbors} = empty
        GrB_TRY (GrB_assign (candidates, new_neighbors, NULL, empty,
            GrB_ALL, n, GrB_DESC_S)) ;

        // sparsify the random number seeds (just keep it for each candidate) 
        // Seed{candidates,replace} = Seed
        GrB_TRY (GrB_assign (Seed, candidates, NULL, Seed, GrB_ALL, n,
            GrB_DESC_RS)) ;

        // Check for stall (can only occur if the matrix has self-edges, or in
        // the exceedingly rare case that 2 nodes have the exact same score).
        // If the method happens to stall, with no nodes selected because
        // the scores happen to tie, try again with another random score.
        GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
        if (last_ncandidates == ncandidates)
        {
            // This case is nearly untestable since it can almost never occur.
            nstall++ ;
            // terminate if the method has stalled too many times
            LG_CHECK (nstall > 32, -111, "stall") ;
            // recreate the random number seeds with a new starting seed
            LAGraph_TRY (LAGraph_Random_Seed (Seed, seed + nstall, msg)) ;
        }
        last_ncandidates = ncandidates ;

        // get the next random Seed vector
        LAGraph_TRY (LAGraph_Random_Next (Seed, msg)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_TRY (LAGraph_Vector_wait (iset, msg)) ;
    (*mis) = iset ;
    LAGraph_FREE_WORK ;
    return (0) ;
}

