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
    GrB_free (&prob) ;              \
    GrB_free (&neighbor_max) ;      \
    GrB_free (&new_members) ;       \
    GrB_free (&new_neighbors) ;     \
    GrB_free (&candidates) ;        \
    GrB_free (&degree) ;            \
    GrB_free (&Seed) ;              \
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
// (diagonal entries) will cause the method to stall.

// Singletons require special treatment.  Since they have no neighbors, their
// prob is never greater than the max of their neighbors, so they never get
// selected and cause the method to stall.  To avoid this case they are removed
// from the candidate set at the begining, and added to the independent set.

// TODO: add an input vector non_candidates, where non_candidates(i) = true
// if node i is not to be considered as a candidate, and is also not added to
// the independent set.

int LAGraph_MaximalIndependentSet       // maximal independent set
(
    // outputs:
    GrB_Vector *mis,            // mis(i) = true if i is in the set
    // inputs:
    LAGraph_Graph G,            // input graph
    int64_t seed,               // random number seed
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector iset = NULL ;
    GrB_Vector prob = NULL ;            // random probability for each node
    GrB_Vector neighbor_max = NULL ;    // value of max neighbor probability
    GrB_Vector new_members = NULL ;     // set of new members to iset
    GrB_Vector new_neighbors = NULL ;   // new neighbors to new iset members
    GrB_Vector candidates = NULL ;      // candidate members to iset
    GrB_Vector empty = NULL ;
    GrB_Vector Seed = NULL ;
    GrB_Vector degree = NULL ;
    GrB_Matrix A ;
    GrB_Index n ;

    LG_CHECK (LAGraph_CheckGraph (G, msg), -102, "graph is invalid") ;
    LG_CHECK (mis == NULL, -103, "mis is null") ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_pattern_is_symmetric == LAGRAPH_TRUE))
    {
        // the pattern of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -105, "G->A must be symmetric") ;
    }

    LG_CHECK (G->rowdegree == NULL, -106, "G->rowdegree must be defined") ;
    LG_CHECK (G->ndiag != 0, -107, "G->ndiag must be zero") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Vector_new (&prob, GrB_FP64, n)) ;
    GrB_TRY (GrB_Vector_new (&neighbor_max, GrB_FP64, n)) ;
    GrB_TRY (GrB_Vector_new (&new_members, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&new_neighbors, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&candidates, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&empty, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&iset, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&Seed, GrB_INT64, n)) ;

    #if LG_SUITESPARSE
    GrB_Semiring symbolic = GxB_ANY_PAIR_BOOL ;
    #else
    GrB_Semiring symbolic = GrB_LOR_LAND_SEMIRING_BOOL ; 
    #endif

    // create the random number seeds
    GrB_TRY (GrB_assign (Seed, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    LAGraph_TRY (LAGraph_Random_Seed (Seed, seed, msg)) ;

    // compute the degree of each nodes in double
    GrB_TRY (GrB_Vector_new (&degree, GrB_FP64, n)) ;
    GrB_TRY (GrB_assign (degree, NULL, NULL, G->rowdegree, GrB_ALL, n, NULL)) ;

    //--------------------------------------------------------------------------
    // remove singletons (nodes of degree zero)
    //--------------------------------------------------------------------------

    GrB_Index nonsingletons = 0 ;
    GrB_TRY (GrB_Vector_nvals (&nonsingletons, degree)) ;
    if (nonsingletons == n)
    { 
        // printf ("no singletons present\n") ;
        // all nodes have degree 1 or more; all nodes are candidates
        GrB_TRY (GrB_assign (candidates, NULL, NULL, (bool) true, GrB_ALL, n,
            NULL)) ;
    }
    else
    {
        // printf ("singletons present: %ld\n", n - nonsingletons) ;
        // one or more singletons are present.  singletons are not candidates;
        // they are added to iset first instead.
        // candidates{degree} = 1
        GrB_TRY (GrB_assign (candidates, degree, NULL, (bool) true, GrB_ALL, n,
            GrB_DESC_S)) ; 
        // add all singletons to iset
        // iset{!degree,replace} = 1
        GrB_TRY (GrB_assign (iset, degree, NULL, (bool) true, GrB_ALL, n,
            GrB_DESC_RSC)) ;
    }

    //--------------------------------------------------------------------------
    // iterate while there are candidates to check
    //--------------------------------------------------------------------------

    GrB_Index ncandidates ;
    GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
    int64_t last_ncandidates = ncandidates ;
    GrB_Index n1 = (GrB_Index) (0.04 * (double) n) ;
    GrB_Index n2 = (GrB_Index) (0.10 * (double) n) ;

    while (ncandidates > 0)
    {
        // printf ("\n ========================================\n") ;
        // GxB_print (candidates, 2) ;

        // sparsify the random number seeds (just keep it for each candidate) 
        // Seed{candidates,replace} = Seed
        GrB_TRY (GrB_assign (Seed, candidates, NULL, Seed, GrB_ALL, n,
            GrB_DESC_RS)) ;
        // GxB_print (Seed, 3) ;

        // prob = random vector with sparsity pattern the same as candidates
        LAGraph_TRY (LAGraph_Random_FP64 (prob, Seed, msg)) ;

        // prob = prob / degree
        GrB_TRY (GrB_eWiseMult (prob, NULL, NULL, GrB_DIV_FP64, prob, degree,
            NULL)) ;
        // GxB_print (prob, 3) ;

        // compute the max probability of all neighbors
        // neighbor_max{candidates,replace} = prob * A
        GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
        if (ncandidates < n1)
        {
            // push
            GrB_TRY (GrB_vxm (neighbor_max, candidates, NULL,
                GrB_MAX_FIRST_SEMIRING_FP64, prob, A, GrB_DESC_RS)) ;
        }
        else
        {
            // pull
            GrB_TRY (GrB_mxv (neighbor_max, candidates, NULL,
                GrB_MAX_SECOND_SEMIRING_FP64, A, prob, GrB_DESC_RS)) ;
        }
        // GxB_print (neighbor_max, 3) ;

        // select node if its probability is > than all its active neighbors
        // new_members = (prob > neighbor_max) using set union so that nodes
        // with no neighbors fall through to the output, as true.
        GrB_TRY (GrB_eWiseAdd (new_members, NULL, NULL, GrB_GT_FP64,
            prob, neighbor_max, NULL)) ;

        // drop explicit zeros from new_members
        #if LG_SUITESPARSE
        #if GxB_IMPLEMENATION >= GxB_VERSION (5,2,0)
        GrB_TRY (GrB_select (new_members, NULL, NULL, GrB_VALUEEQ_BOOL,
            new_members, (bool) true, NULL)) ;
        #else
        GrB_TRY (GxB_select (new_members, NULL, NULL, GxB_NONZERO,
            new_members, NULL, NULL)) ;
        #endif
        #else
        GrB_TRY (GrB_assign (new_members, new_members, NULL, new_members,
            GrB_ALL, n, GrB_DESC_R)) ;
        #endif
        // GxB_print (new_members, 3) ;

        // add new members to independent set
        // iset{new_members} = true
        GrB_TRY (GrB_assign (iset, new_members, NULL, (bool) true,
            GrB_ALL, n, GrB_DESC_S)) ;
        // GxB_print (iset, 3) ;

        // remove new members from set of candidates
        // candidates{new_members} = empty
        GrB_TRY (GrB_assign (candidates, new_members, NULL, empty,
            GrB_ALL, n, GrB_DESC_S)) ;
        // GxB_print (candidates, 3) ;

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
            GrB_TRY (GrB_vxm (new_neighbors, candidates, NULL, symbolic, 
                new_members, A, GrB_DESC_RS)) ;
        }
        else
        {
            // pull
            GrB_TRY (GrB_mxv (new_neighbors, candidates, NULL, symbolic, 
                A, new_members, GrB_DESC_RS)) ;
        }
        // GxB_print (new_neighbors, 3) ;

        // remove new neighbors of new members from set of candidates
        // candidates{new_neighbors} = empty
        GrB_TRY (GrB_assign (candidates, new_neighbors, NULL, empty,
            GrB_ALL, n, GrB_DESC_S)) ;
        // GxB_print (candidates, 3) ;

        GrB_TRY (GrB_Vector_nvals (&ncandidates, candidates)) ;
        LG_CHECK (last_ncandidates == ncandidates, -111, "stall") ;
        last_ncandidates = ncandidates ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    #if GxB_IMPLEMENATION < GxB_VERSION (6,0,0)
    GrB_TRY (GrB_wait (&iset)) ;
    #else
    GrB_TRY (GrB_wait (iset, GrB_MATERIALIZE)) ;
    #endif
    #endif

    (*mis) = iset ;
    LAGraph_FREE_WORK ;
    return (0) ;
}

