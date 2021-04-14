//------------------------------------------------------------------------------
// LAGraph_Property_AT: construct G->AT for a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_ALL GrB_free (&AT) ;

#include "LG_internal.h"

int LAGraph_Property_AT     // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to compute G->AT for
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix AT = NULL ;
    LG_CHECK_INIT (G, msg) ;
    GrB_Matrix A = G->A ;
    LAGraph_Kind kind = G->kind ;

    if (G->AT != NULL || kind == LAGRAPH_ADJACENCY_UNDIRECTED)
    {
        // G->AT already computed, or not needed since A is symmetric
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // G->AT = (G->A)'
    //--------------------------------------------------------------------------

    GrB_Type type ;
    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_TRY (GrB_Matrix_new (&AT, G->A_type, ncols, nrows)) ;
    GrB_TRY (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
    G->AT = AT ;
    G->AT_type = G->A_type;

    return (0) ;
}
