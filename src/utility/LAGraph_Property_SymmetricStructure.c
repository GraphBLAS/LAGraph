//------------------------------------------------------------------------------
// LAGraph_Property_SymmetricStructure: determine G->structure_is_symmetric
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Also computes G->AT if not already computed, if G is not an undirected
// graph and G->A is square.

#define LG_FREE_WORK        \
{                           \
    GrB_free (&S1) ;        \
    GrB_free (&S2) ;        \
    GrB_free (&C) ;         \
}

#include "LG_internal.h"

int LAGraph_Property_SymmetricStructure
(
    // input/output:
    LAGraph_Graph G,    // graph to determine the symmetry of structure of A
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix C = NULL, S1 = NULL, S2 = NULL ;
    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;

    LAGraph_Kind kind = G->kind ;
    if (kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        // assume A is symmetric for an undirected graph
        G->structure_is_symmetric = true ;
        return (GrB_SUCCESS) ;
    }

    if (G->structure_is_symmetric != LAGRAPH_UNKNOWN)
    {
        // symmetric property is already known
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // determine the size of A
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Index n, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    if (n != ncols)
    {
        // A is rectangular and thus cannot be symmetric
        G->structure_is_symmetric = false ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // compute the transpose, if not already computed
    //--------------------------------------------------------------------------

    if (G->AT == NULL)
    {
        LG_TRY (LAGraph_Property_AT (G, msg)) ;
    }

    //--------------------------------------------------------------------------
    // check if the structure of A and AT are the same
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;

    // C(i,j) = 1 if both A(i,j) and AT(i,j) exist
    GRB_TRY (GrB_eWiseMult (C, NULL, NULL, GrB_ONEB_BOOL, A, G->AT, NULL)) ;

    GrB_Index nvals1, nvals2 ;
    GRB_TRY (GrB_Matrix_nvals (&nvals1, C)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals2, A)) ;
    G->structure_is_symmetric = (nvals1 == nvals2) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
