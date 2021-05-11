//------------------------------------------------------------------------------
// LAGraph_Property_NDiag: count the # of diagonal entries of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&M) ;         \
    GrB_free (&D) ;         \
    GrB_free (&d) ;         \
}

#include "LG_internal.h"

int LAGraph_Property_NDiag  // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to compute G->ndiag for
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    GrB_Matrix D = NULL ;
    GrB_Matrix M = NULL ;
    GrB_Vector d = NULL ;
    LG_CHECK_INIT (G, msg) ;
    G->ndiag = LAGRAPH_UNKNOWN ;

    //--------------------------------------------------------------------------
    // extract the diagonal and count its entries
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = LAGraph_MIN (nrows, ncols) ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )

        #if ( GxB_IMPLEMENTATION >= GxB_VERSION (5,0,2) )

            //------------------------------------------------------------------
            // SuiteSparse:GraphBLAS v5.0.2: use GxB_Vector_diag
            //------------------------------------------------------------------

            GrB_TRY (GrB_Vector_new (&d, G->A_type, n)) ;
            GrB_TRY (GxB_Vector_diag (d, A, 0, NULL)) ;
            GrB_TRY (GrB_Vector_nvals (&(G->ndiag), d)) ;

        #else

            //------------------------------------------------------------------
            // SuiteSparse:GraphBLAS v4.x: use GxB_Select
            //------------------------------------------------------------------

            GrB_TRY (GrB_Matrix_new (&D, G->A_type, nrows, ncols)) ;
            GrB_TRY (GxB_select (D, NULL, NULL, GxB_DIAG, A, NULL, NULL)) ;
            GrB_TRY (GrB_Matrix_nvals (&(G->ndiag), D)) ;

        #endif

    #else

        //----------------------------------------------------------------------
        // pure GrB version with no GxB extensions
        //----------------------------------------------------------------------

        GrB_TRY (GrB_Matrix_new (&M, GrB_BOOL, nrows, ncols)) ;
        GrB_TRY (GrB_Matrix_new (&D, G->A_type, nrows, ncols)) ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            // M (i,i) = true
            GrB_TRY (GrB_Matrix_setElement (M, (bool) true, i, i)) ;
        }

        // D<M,struct> = A
        GrB_TRY (GrB_assign (D, M, NULL, A, GrB_ALL, nrows, GrB_ALL, ncols,
            GrB_DESC_S)) ;

        GxB_print (D, 3) ;
        GrB_TRY (GrB_Matrix_nvals (&(G->ndiag), D)) ;

    #endif

    LAGraph_FREE_ALL ;
    return (0) ;
}

