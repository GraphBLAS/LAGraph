//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/dot3_demo.c: test GrB_mxm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
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
    GrB_free (&C) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix C = NULL, R = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = true ;               // set true for diagnostic outputs
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
        false,      // make the graph undirected, if true
        false,      // remove self-edges, if true
        false,      // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

//  printf ("\n==========================The input graph matrix G:\n") ;
//  LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // call GrB_mxm many times
    //--------------------------------------------------------------------------

    uint64_t n ;
    GrB_Matrix A = G->A ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new (&C, GrB_UINT32, n, n)) ;
    GRB_TRY (GrB_Matrix_new (&R, GrB_UINT32, n, n)) ;

    // R<A,struct> = A*A'
    GRB_TRY (GrB_mxm (R, A, NULL, LAGraph_plus_one_uint32, A, A,
        GrB_DESC_RST1)) ;

    #define NTRIALS 10000
    printf ("GrB_mxm: C<A> = A*A', %d trials\n", NTRIALS) ;

    t = LAGraph_WallClockTime ( ) ;
    double t1 = t ;
    for (int k = 0 ; k < NTRIALS ; k++)
    {
        // C<A,struct> = A*A'
        GRB_TRY (GrB_mxm (C, A, NULL, LAGraph_plus_one_uint32, A, A,
            GrB_DESC_RST1)) ;
        double tt = LAGraph_WallClockTime ( ) - t1 ;
        #if 0
        if (tt > 3)
        {
            printf ("%d ok, %g sec\n", k, LAGraph_WallClockTime ( ) - t) ;
            fflush (stdout) ;
            t1 = LAGraph_WallClockTime ( ) ;
        }
        #endif
        #ifndef GRAPHBLAS_HAS_CUDA
        // check the result
        GRB_TRY (GxB_Matrix_fprint (C, "C", 0, stdout)) ;
        bool ok ;
        LG_TRY (LAGraph_Matrix_IsEqual (&ok, C, R, msg)) ;
        if (!ok)
        {
            printf ("failure at trial %d\n", k) ;
            fflush (stdout) ;
            abort ( ) ;
        }
        #endif
    }
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for %d trials:          %g sec\n", NTRIALS, t) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

