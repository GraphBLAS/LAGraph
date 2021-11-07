//------------------------------------------------------------------------------
// LAGraph_Property_ASymmetricStructure: determine G->A_structure_is_symmetric
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

int LAGraph_Property_ASymmetricStructure  // 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to determine the symmetry of structure of A
    char *msg
)
{
    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix C = NULL, S1 = NULL, S2 = NULL ;
    LG_CHECK_INIT (G, msg) ;

    LAGraph_Kind kind = G->kind ;
    if (kind == LAGRAPH_ADJACENCY_UNDIRECTED)
    {
        // assume A is symmetric for an undirected graph
        G->A_structure_is_symmetric = true ;
        return (0) ;
    }

    if (G->A_structure_is_symmetric != LAGRAPH_UNKNOWN)
    {
        // symmetric property is already known
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
        G->A_structure_is_symmetric = false ;
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
    // check if the structure of A and AT are the same
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;

    // C(i,j) = 1 if both A(i,j) and AT(i,j) exist
    GrB_TRY (GrB_eWiseMult (C, NULL, NULL, GrB_ONEB_BOOL, A, G->AT, NULL)) ;

    GrB_Index nvals1, nvals2 ;
    GrB_TRY (GrB_Matrix_nvals (&nvals1, C)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals2, A)) ;
    G->A_structure_is_symmetric = (nvals1 == nvals2) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    return (0) ;
}
