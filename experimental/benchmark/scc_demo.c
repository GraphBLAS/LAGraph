//------------------------------------------------------------------------------
// LAGraph/src/benchmark/cc_demo.c: benchmark LAGr_ConnectedComponents
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Usage: test_cc can be used with both stdin or a file as its input,
// in either grb or mtx format.

//------------------------------------------------------------------------------

#include "LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_alg_internal.h"

#undef  LG_FREE_ALL
#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&components) ;        \
    GrB_free (&components2) ;       \
}

// to run just once, with p = omp_get_max_threads() threads
#define NTHREAD_LIST 1
#define THREAD_LIST 0

// to run with p and p/2 threads, if p = omp_get_max_threads()
// #define NTHREAD_LIST 2
// #define THREAD_LIST 0

// #define NTHREAD_LIST 7
// #define THREAD_LIST 32, 24, 16, 8, 4, 2, 1

// #define NTHREAD_LIST 4
// #define THREAD_LIST 32, 24, 16, 8

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

GrB_Index countCC (GrB_Vector f, GrB_Index n)
{
    GrB_Index nCC = 0;
    GrB_Index *w_val = NULL ;
    LAGraph_Malloc ((void **) &w_val, n, sizeof (GrB_Index), NULL) ;
    if (w_val == NULL) { printf ("out of memory\n") ; abort ( ) ; }
    GrB_Index *i_val = NULL ;
    #if LAGRAPH_SUITESPARSE
    // SuiteSparse:GraphBLAS allows NULL inputs to GrB_Vector_extractTuples
    #else
    LAGraph_Malloc ((void **) &i_val, n, sizeof (GrB_Index), NULL) ;
    if (i_val == NULL) { printf ("out of memory\n") ; abort ( ) ; }
    #endif
    GrB_Vector_extractTuples (i_val, w_val, &n, f) ;
    for (GrB_Index i = 0; i < n; i++)
    {
        if (w_val[i] == i)
        {
            nCC++ ;
        }
    }
    LAGraph_Free ((void **) &i_val, NULL) ;
    LAGraph_Free ((void **) &w_val, NULL) ;
    return nCC;
}

int main (int argc, char **argv)
{

    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Vector components = NULL, components2 = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max, nthreads_outer, nthreads_inner ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
    nthreads_max = nthreads_outer * nthreads_inner ;
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
    printf ("\n%s:\n", matrix_name) ;
    LAGRAPH_TRY (readproblem (&G,
        NULL,   // no source nodes
        false,  // make the graph directed
        false,  // do not remove self-edges
        true,   // structural only, no values needed
        NULL,   // no type preference
        false,  // do not ensure all entries positive
        argc, argv)) ;
    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
    fflush (stdout) ; fflush (stderr) ;

    //--------------------------------------------------------------------------
    // begin tests
    //--------------------------------------------------------------------------

    // warmup
    LAGRAPH_TRY (LAGraph_scc (&components, G->A, msg)) ;
    GrB_Index nCC = countCC (components, n) ;
    printf ("nCC: %llu\n", (long long unsigned int) nCC) ;

#if 0 & LG_CHECK_RESULT
    double tcheck = LAGraph_WallClockTime ( ) ;
    int result = LG_check_cc (components, G, msg) ;
    if (result != 0)
    {
        printf ("test failure: (%d) %s\n", result, msg) ;
    }
    tcheck = LAGraph_WallClockTime ( ) - tcheck ;
    LAGRAPH_TRY (result) ;
    printf ("LG_check_cc passed, time: %g\n", tcheck) ;
#endif

    #define NTRIALS 16
    // #define NTRIALS 1
    printf ("# of trials: %d\n\n", NTRIALS) ;
    fflush (stdout) ; fflush (stderr) ;

    //--------------------------------------------------------------------------
    // LAGr_ConnectedComponents
    //--------------------------------------------------------------------------

    for (int trial = 1 ; trial <= nt ; trial++)
    {
        int nthreads = Nthreads [trial] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, NULL)) ;
        double ttt = 0 ;
        int ntrials = NTRIALS ;
        for (int k = 0 ; k < ntrials ; k++)
        {
            GrB_free (&components2) ;
            double ttrial = LAGraph_WallClockTime ( ) ;
            LAGRAPH_TRY (LAGraph_scc (&components2, G->A, msg)) ;
            ttrial = LAGraph_WallClockTime ( ) - ttrial ;
            ttt += ttrial ;
            printf ("SCC:      nthreads: %2d trial: %2d time: %10.4f sec\n",
                nthreads, k, ttrial) ;
            GrB_Index nCC2 = countCC (components2, n) ;
            if (nCC != nCC2) printf ("failure! %g %g diff %g\n",
                (double) nCC, (double) nCC2, (double) (nCC-nCC2)) ;
            fflush (stdout) ; fflush (stderr) ;
        }
        ttt = ttt / ntrials ;

        printf (         "Avg: CC threads %3d: %10.3f sec, graph: %s\n",
                nthreads, ttt, matrix_name) ;
        fprintf (stderr, "Avg: CC threads %3d: %10.3f sec, graph: %s\n",
                nthreads, ttt, matrix_name) ;
        fflush (stdout) ; fflush (stderr) ;

    }
    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
