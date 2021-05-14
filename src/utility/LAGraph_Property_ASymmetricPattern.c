//------------------------------------------------------------------------------
// LAGraph_Property_ASymmetricPattern: determine G->A_pattern_is_symmetric
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// Also computes G->AT if not already computed, if G is not an undirected
// graph and G->A is square.

#define LAGraph_FREE_WORK   \
{                           \
    GrB_free (&S1) ;        \
    GrB_free (&S2) ;        \
    GrB_free (&C) ;         \
}

#include "LG_internal.h"

int LAGraph_Property_ASymmetricPattern  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine the symmetry of pattern of A
    char *msg
)
{
    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix C = NULL, S1 = NULL, S2 = NULL ;
    LG_CHECK_INIT (G, msg) ;

    G->A_pattern_is_symmetric = LAGRAPH_UNKNOWN ;
    LAGraph_Kind kind = G->kind ;
    if (kind == LAGRAPH_ADJACENCY_UNDIRECTED)
    {
        // assume A is symmetric for an undirected graph
        G->A_pattern_is_symmetric = true ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // determine the size of A
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Index n, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    if (n != ncols)
    {
        // A is rectangular and thus cannot be symmetric
        G->A_pattern_is_symmetric = false ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // compute the transpose, if not already computed
    //--------------------------------------------------------------------------

    if (G->AT == NULL)
    {
        LAGraph_TRY (LAGraph_Property_AT (G, msg)) ;
    }

    //--------------------------------------------------------------------------
    // check if the pattern of A and AT are the same
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;

    #if LG_SUITESPARSE

        // C(i,j) = 1 if both A(i,j) and AT(i,j) exist
        GrB_TRY (GrB_eWiseMult (C, NULL, NULL, GxB_PAIR_BOOL, A, G->AT, NULL)) ;

    #else

        // S1 = pattern of A, with S1(i,j)=1 if A(i,j) exists
        LAGraph_TRY (LAGraph_Pattern (&S1, A, msg)) ;
        // S2 = pattern of AT, with S2(i,j)=1 if AT(i,j) exists
        LAGraph_TRY (LAGraph_Pattern (&S2, G->AT, msg)) ;
        // C(i,j) = 1 if both S1(i,j) == 1 and S2(i,j) == 1
        GrB_TRY (GrB_eWiseMult (C, NULL, NULL, GrB_LAND, S1, S2, NULL)) ;

    #endif

    GrB_Index nvals1, nvals2 ;
    GrB_TRY (GrB_Matrix_nvals (&nvals1, C)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals2, A)) ;
    G->A_pattern_is_symmetric = (nvals1 == nvals2) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    return (0) ;
}
