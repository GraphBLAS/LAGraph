//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/jaccard_demo.c: a demo of the Jaccard algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Elaheh Hassani and Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// This main program is a simple driver for testing and benchmarking the
// LAGr_Jaccard algorithm, in experimental/algorithm.  To use it,
// compile LAGraph while in the build folder with these commands:
//
//      cd LAGraph/build
//      cmake ..
//      make -j8
//
// Then run this demo with an input matrix.  For example:
//
//      ./experimental/benchmark/jaccard_demo ../data/west0067.mtx
//      ./experimental/benchmark/jaccard_demo < ../data/west0067.mtx
//      ./experimental/benchmark/jaccard_demo ../data/karate.mtx
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
    GrB_Matrix JC = NULL ;

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
        true,       // make the graph undirected, if true
        true,       // remove self-edges, if true
        true,       // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program

    t = LAGraph_WallClockTime ( ) - t ;

    int outer, inner ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&outer, &inner, msg)) ;
    printf ("threads (default): %d, %d\n", outer, inner) ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LG_TRY(LAGraph_Graph_Print (G, 1, stdout, msg));
    LG_TRY(LAGraph_Cached_OutDegree(G, msg));

    burble = true ;
    LG_SET_BURBLE (burble) ;

    for (int all_pairs = 0 ; all_pairs <= 1 ; all_pairs++)
    {
        for (int nthreads = inner ; nthreads >= 1 ; )
        {
            if (burble)
            {
                printf ("\n--------------- nthreads %d, all_pairs %d ---------\n",
                    nthreads, all_pairs) ;
            }
            LAGRAPH_TRY (LAGraph_SetNumThreads (outer, nthreads, msg)) ;
            t = LAGraph_WallClockTime ( ) ;
            LG_TRY (LAGr_Jaccard (&JC, G, (bool) all_pairs, msg)) ;
            t = LAGraph_WallClockTime ( ) - t ;
            GrB_free (&JC) ;
            printf ("Time for LAGr_Jaccard (all_pairs %d), nthreads %2d: %g sec\n",
                all_pairs, nthreads, t) ;
            if (nthreads == 32)
            {
                nthreads = 24 ;
            }
            else if (nthreads == 24)
            {
                nthreads = 16 ;
            }
            else
            {
                nthreads = nthreads / 2 ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

