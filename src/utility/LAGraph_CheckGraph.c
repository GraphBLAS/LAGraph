//------------------------------------------------------------------------------
// LAGraph_CheckGraph: check if a graph is valid
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_CheckGraph      // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to check
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear the msg and check basic components
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG_AND_BASIC_ASSERT (G, msg) ;
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
        LG_ASSERT_MSG (nrows == ncols, -1101,
            "adjacency matrix must be square") ;
    }

    #if LG_SUITESPARSE
        // only by-row format is supported when using SuiteSparse
        GxB_Format_Value fmt ;
        GrB_TRY (GxB_get (A, GxB_FORMAT, &fmt)) ;
        LG_ASSERT_MSG (fmt == GxB_BY_ROW, -1102,
            "only by-row format supported") ;
    #endif

    //--------------------------------------------------------------------------
    // check the cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix AT = G->AT ;
    if (AT != NULL)
    {
        GrB_Index nrows2, ncols2;
        GrB_TRY (GrB_Matrix_nrows (&nrows2, AT)) ;
        GrB_TRY (GrB_Matrix_ncols (&ncols2, AT)) ;
        LG_ASSERT_MSG (nrows == ncols2 && ncols == nrows2, -1103,
            "G->AT matrix has the wrong dimensions") ;

        #if LG_SUITESPARSE
            // only by-row format is supported when using SuiteSparse
            GxB_Format_Value fmt ;
            GrB_TRY (GxB_get (AT, GxB_FORMAT, &fmt)) ;
            LG_ASSERT_MSG (fmt == GxB_BY_ROW, -1104,
                "only by-row format supported") ;
        #endif

        // ensure the types of A and AT are the same
        char atype [LAGRAPH_MAX_NAME_LEN] ;
        char ttype [LAGRAPH_MAX_NAME_LEN] ;
        LG_TRY (LAGraph_MatrixTypeName (atype, A, msg)) ;
        LG_TRY (LAGraph_MatrixTypeName (ttype, AT, msg)) ;
        LG_ASSERT_MSG (MATCHNAME (atype, ttype),
            -1105, "A and AT must have the same type") ;
    }

    GrB_Vector rowdegree = G->rowdegree ;
    if (rowdegree != NULL)
    {
        GrB_Index m ;
        GrB_TRY (GrB_Vector_size (&m, rowdegree)) ;
        LG_ASSERT_MSG (m == nrows, -1106, "rowdegree invalid size") ;
        char rtype [LAGRAPH_MAX_NAME_LEN] ;
        LG_TRY (LAGraph_VectorTypeName (rtype, rowdegree, msg)) ;
        LG_ASSERT_MSG (MATCHNAME (rtype, "int64_t"),
            -1107, "rowdegree has wrong type; must be GrB_INT64") ;
    }

    GrB_Vector coldegree = G->coldegree ;
    if (coldegree != NULL)
    {
        GrB_Index n ;
        GrB_TRY (GrB_Vector_size (&n, coldegree)) ;
        LG_ASSERT_MSG (n == ncols, -1108, "coldegree invalid size") ;
        char ctype [LAGRAPH_MAX_NAME_LEN] ;
        LG_TRY (LAGraph_VectorTypeName (ctype, coldegree, msg)) ;
        LG_ASSERT_MSG (MATCHNAME (ctype, "int64_t"),
            -1109, "coldegree has wrong type; must be GrB_INT64") ;
    }

    return (GrB_SUCCESS) ;
}
