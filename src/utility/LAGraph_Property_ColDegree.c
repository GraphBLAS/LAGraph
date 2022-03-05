//------------------------------------------------------------------------------
// LAGraph_Property_ColDegree: determine G->coldegree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_Property_ColDegree computes G->coldegree, where G->coldegree(j) is
// the number of entries in G->A (:,j).  If there are no entries in G->A (:,j),
// G->coldgree(j) is not present in the structure of G->coldegree.  That is,
// G->coldegree contains no explicit zero entries.

// G->coldegree is not computed if the graph is undirected.  Use G->rowdegree
// instead, and LAGraph_Property_RowDegree.

#define LG_FREE_WORK            \
{                               \
    GrB_free (&S) ;             \
    GrB_free (&x) ;             \
}

#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
    GrB_free (&coldegree) ;     \
}

#include "LG_internal.h"

int LAGraph_Property_ColDegree
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->coldegree
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix S = NULL ;
    GrB_Vector coldegree = NULL, x = NULL ;
    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;

    if (G->coldegree != NULL || G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        // G->coldegree already computed, or not needed
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // determine the size of A
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Matrix AT = G->AT ;
    GrB_Index nrows, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

    //--------------------------------------------------------------------------
    // compute the coldegree
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&coldegree, GrB_INT64, ncols)) ;
    // x = zeros (nrows,1)
    GRB_TRY (GrB_Vector_new (&x, GrB_INT64, nrows)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, nrows, NULL)) ;

    if (AT != NULL)
    {
        // G->coldegree = row degree of AT; this will be faster assuming
        // AT is held in a row-oriented format. 
        GRB_TRY (GrB_mxv (coldegree, NULL, NULL, LAGraph_plus_one_int64,
            AT, x, NULL)) ;
    }
    else
    {
        // G->coldegree = column degree of A
        GRB_TRY (GrB_mxv (coldegree, NULL, NULL, LAGraph_plus_one_int64,
            A, x, GrB_DESC_T0)) ;
    }

    G->coldegree = coldegree ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
