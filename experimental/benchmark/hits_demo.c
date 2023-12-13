//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/hits_demo.c: a demo of the HITS algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Aurko Routh, Texas A&M University

//------------------------------------------------------------------------------

// This main program is a simple driver for testing and benchmarking the
// LAGr_HITS "algorithm", in experimental/algorithm.  To use it,
// compile LAGraph while in the build folder with these commands:
//
//      cd LAGraph/build
//      cmake ..
//      make -j8
//
// Then run this demo with an input matrix.  For example:
//
//      ./experimental/benchmark/hits_demo ../data/west0067.mtx
//      ./experimental/benchmark/hits_demo < ../data/west0067.mtx
//      ./experimental/benchmark/hits_demo ../data/karate.mtx
//
#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Vector hubs = NULL ;
    GrB_Vector authorities = NULL;

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
        false,      // make the graph undirected, if true
        false,      // remove self-edges, if true
        true,      // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program

    t = LAGraph_WallClockTime ( ) - t ;
    LAGRAPH_TRY (LAGraph_SetNumThreads (1, 24, msg)) ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LG_TRY(LAGraph_Cached_OutDegree(G, msg));
    LG_TRY(LAGraph_Cached_InDegree(G, msg));
    
    int iters = 0, itermax = 1000;
    float tol = 1e-6;
    t = LAGraph_WallClockTime ( ) ;
    
    LG_TRY (LAGr_HITS (&hubs, &authorities, &iters, G, tol, itermax, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGr_HITS: %g sec\n", t) ;
    //-------------------------------------------------------------------------- 
    // free everyting and finish
    //--------------------------------------------------------------------------
    printf("Num iterations: %d\n", iters);
    
    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
    