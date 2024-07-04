//------------------------------------------------------------------------------
// LAGraph_RegularReachability.c: regular reachability
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Georgiy Belyanin, Semyon Grigoriev, St. Petersburg State
// University.

//------------------------------------------------------------------------------

// For a labeled directed graph the algorithm computes the set of nodes for
// which these conditions are held:
// * The node is reachable by a path from one of the source nodes.
// * The concatenation of the labels over this path's edges satisfies the
//   specified regular restrictions.
//
// The regular restrictions are specified by a non-determitistic finite
// automaton (NFA) over the same graph label alphabet. The algorithm assumes
// the labels are enumerated from 0 to some nl-1. The input graph and the NFA
// are defined by adjacency matrix decomposition. They are represented by arrays
// of graphs G and R both of length nl in which G[i]/R[i] represents the
// adjacency matrix of the i-th label for the graph/NFA correspondingly.
//
// The algorithm is based on the idea of considering the graph as
// non-deterministic finite automaton having the specified set of stating nodes
// as starting states and evaulating two modified breadth-first traversals of
// the graph and the input NFA at the same time considering the matching labels.
//
// The intuition behind this proccess is similar to building a direct product
// of the graph automaton and the input automaton. This construction results
// with an intersection of two regular languages. The first one is defined by
// the set of all words that can be obtained through edge label concatination
// of paths starting in one of the source nodes. And the second one is the set
// of the paths accepted by the NFA thus matching the desired restrictions.
//
// On algorithm step n the relation between the NFA edges and the graph nodes
// is bulit. These conditions are held for the pairs in this relation:
// * The node is reachable by a path of length n from one of the source nodes
//   in the graph.
// * The state is reachable by a path of length n from one of the starting
//   states in the NFA.
// * These paths have the same labels.
// The algorithm accumulates these relations. Then it extracts the nodes that
// are in relation with one of the final states.
//
// Full description is available at:
//   http://arxiv/will/be/here/i/hope
//
// Performance considerations: the best performance is shown when the algorithm
// receives a minimal deterministic finite automaton as an input.

#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&frontier) ;                      \
    GrB_free (&next_frontier) ;                 \
    GrB_free (&symbol_frontier) ;               \
    LAGraph_Free ((void **) &A, NULL) ;         \
    LAGraph_Free ((void **) &AT, NULL) ;        \
    LAGraph_Free ((void **) &B, NULL) ;         \
    LAGraph_Free ((void **) &BT, NULL) ;        \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (reachable) ;          \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_RegularReachability
(
    // output:
    GrB_Vector *reachable,      // reachable(i) = true if node i is reachable
                                // from one of the starting nodes
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacenc ymatrix decompositions
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA 
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA 
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Matrix frontier = NULL ;
    GrB_Matrix symbol_frontier = NULL ;
    GrB_Matrix next_frontier = NULL ;
    GrB_Matrix visited = NULL ;
    GrB_Vector final_reducer = NULL ;

    GrB_Index ng = 0 ;                   // # nodes in the graph
    GrB_Index nr = 0 ;                   // # states odes in the NFA
    GrB_Index states = ns ;              // # pairs in the current
                                         // correspondence beetwen the graph and
                                         // the NFA

    GrB_Index rows = 0 ;                 // utility matrix row count
    GrB_Index cols = 0 ;                 // utility matrix column count
    GrB_Index vals = 0 ;                 // utility matrix value count

    GrB_Matrix *A = NULL ;
    GrB_Matrix *AT = NULL ;
    GrB_Matrix *B = NULL ;
    GrB_Matrix *BT = NULL ;

    LG_ASSERT (reachable != NULL && G != NULL && R != NULL && S != NULL,
        GrB_NULL_POINTER) ;
    (*reachable) = NULL ;

    for (int32_t i = 0; i < nl; i++)
    {
        if (G[i] == NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (G[i], msg)) ;
    }

    for (int32_t i = 0; i < nl; i++)
    {
        if (R[i] == NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (R[i], msg)) ;
    }

    LG_TRY (LAGraph_Malloc ((void **) &A, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &AT, nl, sizeof (GrB_Matrix), msg)) ;

    for (int32_t i = 0; i < nl; i++) 
    {
        if (G[i] == NULL) 
        {
            A[i] = NULL ;
            continue ;
        }

        A[i] = G[i]->A;
        if (G[i]->kind == LAGraph_ADJACENCY_UNDIRECTED ||
            G[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            AT[i] = A[i] ;
        }
        else
        {
            AT[i] = G[i]->AT ;
            LG_ASSERT_MSG (AT[i] != NULL, LAGRAPH_NOT_CACHED,
                "G[i]->AT is required") ;
        }
    }

    LG_TRY (LAGraph_Malloc ((void **) &B, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &BT, nl, sizeof (GrB_Matrix), msg)) ;

    for (int32_t i = 0; i < nl; i++)
    {
        if (R[i] == NULL)
        {
            B[i] = NULL ;
            continue ;
        }

        B[i] = R[i]->A;
        if (R[i]->kind == LAGraph_ADJACENCY_UNDIRECTED ||
            R[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            BT[i] = B[i] ;
        }
        else
        {
            BT[i] = R[i]->AT ;
            LG_ASSERT_MSG (BT[i] != NULL, LAGRAPH_NOT_CACHED,
                "R[i]->AT is required") ;
        }
    }

    for (int i = 0; i < nl; i++) 
    {
        if (A[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&ng, A[i])) ;
        break ;
    }

    for (int i = 0; i < nl; i++) 
    {
        if (B[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&nr, B[i])) ;
        break ;
    }

    // Check all the matricies in graph adjacency matrix decomposition are
    // square and of the same dimensions
    for (int i = 0; i < nl; i++) 
    {
        if (A[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&rows, A[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, A[i])) ;

        LG_ASSERT_MSG (rows == ng && cols == ng, LAGRAPH_NOT_CACHED,
            "all the matricies in the graph adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check all the matricies in NFA adjacency matrix decomposition are
    // square and of the same dimensions
    for (int i = 0; i < nl; i++) 
    {
        if (B[i] == NULL) continue ;

        GrB_Index rows = 0 ;
        GrB_Index cols = 0 ;

        GRB_TRY (GrB_Matrix_nrows (&rows, B[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, B[i])) ;

        LG_ASSERT_MSG (rows == nr && cols == nr, LAGRAPH_NOT_CACHED,
            "all the matricies in the NFA adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check source nodes in the graph
    for (GrB_Index i = 0; i < ns; i++)
    {
        GrB_Index s = S [i] ;
        LG_ASSERT_MSG (s < ng, GrB_INVALID_INDEX, "invalid graph source node") ;
    }

    // Check starting states of the NFA
    for (GrB_Index i = 0; i < nqs; i++)
    {
        GrB_Index qs = QS [i] ;
        LG_ASSERT_MSG (qs < nr, GrB_INVALID_INDEX,
            "invalid NFA starting state") ;
    }

    // Check final states of the NFA
    for (GrB_Index i = 0; i < nqf; i++)
    {
        GrB_Index qf = QF [i] ;
        LG_ASSERT_MSG (qf < nr, GrB_INVALID_INDEX, "invalid NFA final state") ;
    }

    // -------------------------------------------------------------------------
    // initialization
    // -------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (reachable, GrB_BOOL, ng)) ;

    GRB_TRY (GrB_Vector_new (&final_reducer, GrB_BOOL, nr)) ;

    // Initialize matrix for reducing the result
    for (GrB_Index i = 0; i < nqf; i++)
    {
        GrB_Index qf = QF[i] ;

        GRB_TRY (GrB_Vector_setElement (final_reducer, true, qf)) ;
    }

    GRB_TRY (GrB_Matrix_new (&next_frontier, GrB_BOOL, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&visited, GrB_BOOL, nr, ng)) ;

    // Initialize frontier with the source nodes
    for (GrB_Index i = 0; i < nqs; i++)
    {
        GrB_Index qs = QS[i] ;

        for (GrB_Index j = 0; j < ns; j++)
        {
            GrB_Index s = S [j] ;
            GRB_TRY (GrB_Matrix_setElement (next_frontier, true, qs, s)) ;
            GRB_TRY (GrB_Matrix_setElement (visited, true, qs, s)) ;
        }
    }

    // Initialize a few utility matricies
    GRB_TRY (GrB_Matrix_new (&frontier, GrB_BOOL, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&symbol_frontier, GrB_BOOL, nr, ng)) ;

    // Main loop
    while (states != 0)
    {
        GrB_Matrix old_frontier = frontier ;
        frontier = next_frontier ;
        next_frontier = old_frontier ;

        GRB_TRY (GrB_Matrix_clear(next_frontier)) ;

        // Obtain a new relation between the NFA states and the graph nodes
        for (int32_t i = 0; i < nl; i++)
        {
            if (A[i] == NULL || B[i] == NULL) continue ;

            // Traverse the NFA
            GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL,
                GrB_LOR_LAND_SEMIRING_BOOL, BT[i], frontier, GrB_DESC_R)) ;

            // Traverse the graph
            GRB_TRY (GrB_mxm (next_frontier, visited, GrB_LOR,
                GrB_LOR_LAND_SEMIRING_BOOL, symbol_frontier, A[i], GrB_DESC_SC)) ;
        }

        // Accumulate the new state <-> node correspondence
        GRB_TRY (GrB_assign (visited, visited, GrB_NULL, next_frontier,
            GrB_ALL, nr, GrB_ALL, ng, GrB_DESC_SC)) ;

        GRB_TRY (GrB_Matrix_nvals (&states, next_frontier)) ;
    }

    // Extract the nodes matching the final NFA states
    GRB_TRY (GrB_vxm (*reachable, GrB_NULL, GrB_NULL,
        GrB_LOR_LAND_SEMIRING_BOOL, final_reducer, visited, GrB_NULL)) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
