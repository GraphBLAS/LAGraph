//------------------------------------------------------------------------------
// LAGraph_Cached_ColDegree: determine G->col_degree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_Cached_ColDegree computes G->col_degree, where G->col_degree(j) is
// the number of entries in G->A (:,j).  If there are no entries in G->A (:,j),
// G->coldgree(j) is not present in the structure of G->col_degree.  That is,
// G->col_degree contains no explicit zero entries.

// G->col_degree is not computed if the graph is undirected.  Use G->row_degree
// instead, and LAGraph_Cached_RowDegree.

#define LG_FREE_WORK            \
{                               \
    GrB_free (&S) ;             \
    GrB_free (&x) ;             \
}

#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
    GrB_free (&col_degree) ;    \
}

#include "LG_internal.h"

int LAGraph_Cached_ColDegree
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->col_degree
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix S = NULL ;
    GrB_Vector col_degree = NULL, x = NULL ;
    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;

    if (G->col_degree != NULL)
    {
        // G->col_degree already computed
        return (GrB_SUCCESS) ;
    }

    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        // G->col_degree is not computed since A is symmetric (warning only)
        return (LAGRAPH_CACHE_NOT_NEEDED) ;
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
    // compute the col_degree
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&col_degree, GrB_INT64, ncols)) ;
    // x = zeros (nrows,1)
    GRB_TRY (GrB_Vector_new (&x, GrB_INT64, nrows)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, nrows, NULL)) ;

    if (AT != NULL)
    {
        // G->col_degree = row degree of AT; this will be faster assuming
        // AT is held in a row-oriented format. 
        GRB_TRY (GrB_mxv (col_degree, NULL, NULL, LAGraph_plus_one_int64,
            AT, x, NULL)) ;
    }
    else
    {
        // G->col_degree = column degree of A
        GRB_TRY (GrB_mxv (col_degree, NULL, NULL, LAGraph_plus_one_int64,
            A, x, GrB_DESC_T0)) ;
    }

    G->col_degree = col_degree ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
