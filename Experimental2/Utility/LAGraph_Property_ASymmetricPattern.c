//------------------------------------------------------------------------------
// LAGraph_Property_ASymmetricPattern: determine G->A_pattern_is_symmetric
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Also computes G->AT if not already computed, if G is not an undirected
// graph and G->A is square.

#define LAGraph_FREE_WORK GrB_free (&C) ;

#include "LAGraph_Internal.h"

int LAGraph_Property_ASymmetricPattern  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine the symmetry of pattern of A
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix C = NULL ;
    LAGraph_CHECK_INIT (G, msg) ;
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
    GrB_TRY (GrB_eWiseMult (C, NULL, NULL, GxB_PAIR_BOOL, A, G->AT, NULL)) ;

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

