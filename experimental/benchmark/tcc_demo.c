//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/tcc_demo.c: benchmark for triangle centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Tim Davis, Texas A&M

//------------------------------------------------------------------------------

// Usage:  tcc_demo < matrixmarketfile.mtx
//         tcc_demo matrixmarketfile.mtx
//         tcc_demo matrixmarketfile.grb

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

// #define NTHREAD_LIST 2

#define NTHREAD_LIST 1
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

// #define NTHREAD_LIST 7
// #define THREAD_LIST 40, 20, 16, 8, 4, 2, 1

#define LG_FREE_ALL                 \
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
    ntrials = 3 ;
    printf ("# of trials: %d\n", ntrials) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
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
    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, true, NULL, false, argc, argv)) ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
    // LAGRAPH_TRY (LAGraph_DisplayGraph (G, LAGraph_SHORT, stdout, msg)) ;
    // ensure G->A is FP64 and all 1
    GRB_TRY (GrB_Matrix_new (&A, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_assign (A, G->A, NULL, (double) 1,
        GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
    GrB_free (&(G->A)) ;
    G->A = A ;

    //--------------------------------------------------------------------------
    // triangle centrality
    //--------------------------------------------------------------------------

    // warmup for more accurate timing
    double tic [2], tt ;
    uint64_t ntri ;
    LAGRAPH_TRY (LAGraph_Tic (tic, NULL)) ;
    LAGRAPH_TRY (LAGraph_VertexCentrality_Triangle (&c, &ntri, 3, G, msg)) ;
    LAGRAPH_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
    GRB_TRY (GrB_free (&c)) ;
    printf ("warmup time %g sec, # triangles: %g\n", tt, (double) ntri) ;

    for (int method = 1 ; method <= 3 ; method += 2)
    {
        for (int t = 1 ; t <= nt ; t++)
        {
            int nthreads = Nthreads [t] ;
            if (nthreads > nthreads_max) continue ;
            LAGRAPH_TRY (LAGraph_SetNumThreads (nthreads, msg)) ;
            double ttot = 0, ttrial [100] ;
            for (int trial = 0 ; trial < ntrials ; trial++)
            {
                LAGRAPH_TRY (LAGraph_Tic (tic, NULL)) ;
                LAGRAPH_TRY (LAGraph_VertexCentrality_Triangle (&c, &ntri,
                    method, G, msg)) ;
                GRB_TRY (GrB_free (&c)) ;
                LAGRAPH_TRY (LAGraph_Toc (&ttrial [trial], tic, NULL)) ;
                ttot += ttrial [trial] ;
                printf ("threads %2d trial %2d: %12.6f sec\n",
                    nthreads, trial, ttrial [trial]) ;
                fprintf (stderr, "threads %2d trial %2d: %12.6f sec\n", 
                    nthreads, trial, ttrial [trial]) ;
            }
            ttot = ttot / ntrials ;

            printf ("Avg: TCentrality(%d) "
                "nthreads: %3d time: %12.6f matrix: %s\n",
                method, nthreads, ttot, matrix_name) ;

            fprintf (stderr, "Avg: TCentrality(%d) "
                "nthreads: %3d time: %12.6f matrix: %s\n",
                method, nthreads, ttot, matrix_name) ;
        }
    }

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

