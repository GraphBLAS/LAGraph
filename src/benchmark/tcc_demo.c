//------------------------------------------------------------------------------
// LAGraph/src/benchmark/tcc_demo.c: benchmark for triangle centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage:  tcc_demo < matrixmarketfile.mtx
//         tcc_demo matrixmarketfile.mtx
//         tcc_demo matrixmarketfile.grb

#include "LAGraph_demo.h"

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&A) ;                 \
    GrB_free (&c) ;                 \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;

    GrB_Vector c = NULL ;
    GrB_Matrix A = NULL ;
    LAGraph_Graph G = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    int ntrials = 3 ;
    // ntrials = 1 ;        // HACK
    printf ("# of trials: %d\n", ntrials) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ;
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
    if (Nthreads [1] == 0)
    {
        // create thread list automatically
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test: ") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    if (readproblem (&G, NULL,
        true, true, true, NULL, false, argc, argv) != 0) ERROR ;

    GrB_Index n, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // triangle centrality
    //--------------------------------------------------------------------------

    // warmup for more accurate timing
    double tic [2] ;
    LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
    LAGraph_TRY (LAGraph_VertexCentrality_Triangle (&c, G, msg)) ;
    GrB_TRY (GrB_free (&c)) ;

    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_TRY (LAGraph_SetNumThreads (nthreads, msg)) ;
        double ttot = 0, ttrial [100] ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            printf ("\nstart trial %d\n", trial) ;
            LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;

            LAGraph_TRY (LAGraph_VertexCentrality_Triangle (&c, G, msg)) ;
            GrB_TRY (GrB_free (&c)) ;

            LAGraph_TRY (LAGraph_Toc (&ttrial [trial], tic, NULL)) ;
            ttot += ttrial [trial] ;
            printf ("trial %2d: %12.6f sec rate %6.2f\n",
                trial, ttrial [trial], 1e-6 * nvals / ttrial [trial]) ;
        }
        ttot = ttot / ntrials ;

        printf ("nthreads: %3d time: %12.6f rate: %6.2f\n", nthreads,
                ttot, 1e-6 * nvals / ttot) ;
        fprintf (stderr, "Avg: TCentrality threads: %2d  time: %10.3f sec, "
            "matrix: %s\n", nthreads, ttot, matrix_name) ;
    }

    LAGraph_FREE_ALL ;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}

