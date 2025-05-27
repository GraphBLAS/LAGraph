//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/Incidence_Matrix_demo.c: a demo to check the 
// speed of Incidence Matrix building
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&Y) ;                             \
    GrB_free (&E) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix Y = NULL ;
    GrB_Matrix E = NULL ;
    GrB_Type atype = NULL;
    GrB_Semiring op = NULL;
    GrB_Index nrows, ncols;

    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).

    double t = LAGraph_WallClockTime ( ) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        true,      // make the graph undirected, if true
        true,      // remove self-edges, if true
        false,      // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // try the LAGraph_Incidence_Matrix "algorithm"
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    LG_TRY (LAGraph_Incidence_Matrix (&E, G, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_Incidence_Matrix: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results (make sure Y is a copy of G->A)
    //--------------------------------------------------------------------------

    bool isequal ;
    t = LAGraph_WallClockTime ( ) ;
    //unnessesary work but gets me dimensions and type quickly.
    GRB_TRY (GxB_Matrix_type(&atype, G->A)) ;
    GRB_TRY (GrB_Matrix_nrows(&nrows, G->A)) ;
    GRB_TRY (GrB_Matrix_nrows(&ncols, G->A)) ;
    GRB_TRY (GrB_Matrix_new(&Y, atype, nrows, ncols)) ; 
    if      (atype == GrB_BOOL  ) op = GxB_ANY_FIRST_BOOL   ;
    else if (atype == GrB_INT8  ) op = GxB_ANY_FIRST_INT8   ;
    else if (atype == GrB_INT16 ) op = GxB_ANY_FIRST_INT16  ;
    else if (atype == GrB_INT32 ) op = GxB_ANY_FIRST_INT32  ;
    else if (atype == GrB_INT64 ) op = GxB_ANY_FIRST_INT64  ;
    else if (atype == GrB_UINT8 ) op = GxB_ANY_FIRST_UINT8  ;
    else if (atype == GrB_UINT16) op = GxB_ANY_FIRST_UINT16 ;
    else if (atype == GrB_UINT32) op = GxB_ANY_FIRST_UINT32 ;
    else if (atype == GrB_UINT64) op = GxB_ANY_FIRST_UINT64 ;
    else if (atype == GrB_FP32  ) op = GxB_ANY_FIRST_FP32   ;
    else if (atype == GrB_FP64  ) op = GxB_ANY_FIRST_FP64   ;
    GRB_TRY (GrB_mxm (Y, NULL, NULL, op, E, E, GrB_DESC_T1));
    GRB_TRY (GrB_select (Y, NULL, NULL, GrB_OFFDIAG, Y, 0, NULL))
    LG_TRY (LAGraph_Matrix_IsEqual (&isequal, Y, G->A, msg));

    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;
    if (isequal)
    {
        printf ("Test passed.\n") ;
    }
    else
    {
        printf ("Test failure!\n") ;
    }

    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf ("\n===============================The result matrix Y:\n") ;
    LG_TRY (LAGraph_Matrix_Print (Y, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
