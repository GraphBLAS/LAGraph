//------------------------------------------------------------------------------
// LAGraph_CheckGraph: check if a graph is valid
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_Internal.h"

int LAGraph_CheckGraph      // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to check
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear the msg and check basic components
    //--------------------------------------------------------------------------

    LAGraph_CHECK_INIT (G, msg) ;
    GrB_Matrix A = G->A ;
    LAGraph_Kind kind = G->kind ;

    //--------------------------------------------------------------------------
    // ensure the matrix is square for directed or undirected graphs
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    if (kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
        kind == LAGRAPH_ADJACENCY_DIRECTED)
    {
        GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
        GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
        LAGraph_CHECK (nrows != ncols, -1, "adjacency matrix invalid") ;
    }

    // only by-row format is supported
    GxB_Format_Value fmt ;
    GrB_TRY (GxB_get (A, GxB_FORMAT, &fmt)) ;
    LAGraph_CHECK (fmt != GxB_BY_ROW, -1, "only by-row format supported") ;

    //--------------------------------------------------------------------------
    // check the cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix AT = G->AT ;
    if (AT != NULL)
    {
        GrB_Index nrows2, ncols2;
        GrB_TRY (GrB_Matrix_nrows (&nrows2, AT)) ;
        GrB_TRY (GrB_Matrix_ncols (&ncols2, AT)) ;
        LAGraph_CHECK (nrows != ncols2 || ncols != nrows2, -1,
            "G->AT matrix invalid") ;

        // only by-row format is supported
        GxB_Format_Value fmt ;
        GrB_TRY (GxB_get (AT, GxB_FORMAT, &fmt)) ;
        LAGraph_CHECK (fmt != GxB_BY_ROW, -1, "only by-row format supported") ;

        // ensure the types of A and AT are the same
        GrB_Type type1, type2 ;
        GrB_TRY (GxB_Matrix_type (&type1, A)) ;
        GrB_TRY (GxB_Matrix_type (&type2, AT)) ;
        LAGraph_CHECK (type1 != type2, -1, "A and AT types are different") ;
    }

    GrB_Vector rowdegree = G->rowdegree ;
    if (rowdegree != NULL)
    {
        GrB_Index m ;
        GrB_TRY (GrB_Vector_size (&m, rowdegree)) ;
        LAGraph_CHECK (m != nrows, -1, "rowdegree invalid") ;
        GrB_Type type ;
        GrB_TRY (GxB_Vector_type (&type, rowdegree)) ;
        LAGraph_CHECK (type != GrB_INT64, -1, "rowdegree has wrong type") ;
    }

    GrB_Vector coldegree = G->coldegree ;
    if (coldegree != NULL)
    {
        GrB_Index n ;
        GrB_TRY (GrB_Vector_size (&n, coldegree)) ;
        LAGraph_CHECK (n != ncols, -1, "coldegree invalid") ;
        GrB_Type type ;
        GrB_TRY (GxB_Vector_type (&type, coldegree)) ;
        LAGraph_CHECK (type != GrB_INT64, -1, "coldegree has wrong type") ;
    }

    return (0) ;
}

