//------------------------------------------------------------------------------
// LAGraph_Property_RowDegree: determine G->row_degree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_Property_RowDegree computes G->row_degree, where G->row_degree(i) is
// the number of entries in G->A (i,:).  If there are no entries in G->A (i,:),
// G->rowdgree(i) is not present in the structure of G->row_degree.  That is,
// G->row_degree contains no explicit zero entries.

#define LG_FREE_WORK            \
{                               \
    GrB_free (&x) ;             \
}

#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
    GrB_free (&row_degree) ;     \
}

#include "LG_internal.h"

int LAGraph_Property_RowDegree
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->row_degree
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Vector row_degree = NULL, x = NULL ;
    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;

    if (G->row_degree != NULL)
    {
        // G->row_degree already computed
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // determine the size of A
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Index nrows, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

    //--------------------------------------------------------------------------
    // compute the row_degree
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&row_degree, GrB_INT64, nrows)) ;
    // x = zeros (ncols,1)
    GRB_TRY (GrB_Vector_new (&x, GrB_INT64, ncols)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, ncols, NULL)) ;

    GRB_TRY (GrB_mxv (row_degree, NULL, NULL, LAGraph_plus_one_int64,
        A, x, NULL)) ;

    G->row_degree = row_degree ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
