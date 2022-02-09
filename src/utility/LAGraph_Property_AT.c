//------------------------------------------------------------------------------
// LAGraph_Property_AT: construct G->AT for a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LG_FREE_ALL GrB_free (&AT) ;

#include "LG_internal.h"

int LAGraph_Property_AT
(
    // input/output:
    LAGraph_Graph G,    // graph for which to compute G->AT
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix AT = NULL ;
    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;
    GrB_Matrix A = G->A ;
    LAGraph_Kind kind = G->kind ;

    if (G->AT != NULL || kind == LAGRAPH_ADJACENCY_UNDIRECTED)
    {
        // G->AT already computed, or not needed since A is symmetric
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // G->AT = (G->A)'
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Type atype ;
    char atype_name [LAGRAPH_MAX_NAME_LEN] ;
    LG_TRY (LAGraph_MatrixTypeName (atype_name, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&atype, atype_name, msg)) ;
    GrB_TRY (GrB_Matrix_new (&AT, atype, ncols, nrows)) ;
    GrB_TRY (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
    G->AT = AT ;

    return (GrB_SUCCESS) ;
}
