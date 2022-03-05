//------------------------------------------------------------------------------
// LAGraph_Property_RowDegree: determine G->rowdegree
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_Property_RowDegree computes G->rowdegree, where G->rowdegree(i) is
// the number of entries in G->A (i,:).  If there are no entries in G->A (i,:),
// G->rowdgree(i) is not present in the structure of G->rowdegree.  That is,
// G->rowdegree contains no explicit zero entries.

#define LG_FREE_WORK            \
{                               \
    GrB_free (&x) ;             \
}

#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
    GrB_free (&rowdegree) ;     \
}

#include "LG_internal.h"

int LAGraph_Property_RowDegree
(
    // input/output:
    LAGraph_Graph G,    // graph to determine G->rowdegree
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Vector rowdegree = NULL, x = NULL ;
    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;

    if (G->rowdegree != NULL)
    {
        // G->rowdegree already computed
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
    // compute the rowdegree
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&rowdegree, GrB_INT64, nrows)) ;
    // x = zeros (ncols,1)
    GRB_TRY (GrB_Vector_new (&x, GrB_INT64, ncols)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, ncols, NULL)) ;

    GRB_TRY (GrB_mxv (rowdegree, NULL, NULL, LAGraph_plus_one_int64,
        A, x, NULL)) ;

    G->rowdegree = rowdegree ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
