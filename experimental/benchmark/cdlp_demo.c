//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/lcc_demo.c:
// benchmark for community detection using label propagation
//------------------------------------------------------------------------------

// LAGraph, (c) 2023 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Pascal Costanza, Intel, Belgium
// Based on tcc_demo by Tim Davis, Texas A&M

//------------------------------------------------------------------------------

// Usage:  cdlp_demo < matrixmarketfile.mtx
//         cdlp_demo matrixmarketfile.mtx
//         cdlp_demo matrixmarketfile.grb

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_Xtest.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0

#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&c) ;                 \
    GrB_free (&cgood) ;             \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;

    GrB_Vector c = NULL, cgood = NULL ;
    LAGraph_Graph G = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    int ntrials = 3 ;
    ntrials = 3 ;
    printf ("# of trials: %d\n", ntrials) ;

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
    LAGRAPH_TRY (readproblem (&G, NULL,
                              false, true, true, NULL, false, argc, argv)) ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    LAGRAPH_TRY (LAGraph_Cached_IsSymmetricStructure(G, msg)) ;

    //--------------------------------------------------------------------------
    // commmunity detection using label propagation
    //--------------------------------------------------------------------------

    // compute check result
    double tt = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_cdlp_withsort (&cgood, G, 100, msg)) ;
    tt = LAGraph_WallClockTime() - tt ;
    printf ("compute check time %g sec\n", tt) ;

    // warmup for more accurate timing
    tt = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_cdlp (&c, G, 100, msg)) ;
    tt = LAGraph_WallClockTime ( ) - tt ;
    printf ("warmup time %g sec\n", tt) ;

    // check result
    LAGRAPH_TRY (GrB_wait (c, GrB_MATERIALIZE)) ;
    LAGRAPH_TRY (GrB_wait (cgood, GrB_MATERIALIZE)) ;
    bool isequal ;
    LAGRAPH_TRY( LAGraph_Vector_IsEqual (&isequal, c, cgood, msg)) ;
    if (!isequal) {
        LAGRAPH_TRY (GrB_PANIC) ;
    }
    LAGRAPH_TRY (GrB_free (&c)) ;
    LAGRAPH_TRY (GrB_free (&cgood)) ;


    for (int t = 1; t <= nt; t++)
    {
        int nthreads = Nthreads[t];
        if (nthreads > nthreads_max) continue;
        LAGRAPH_TRY (LAGraph_SetNumThreads(1, nthreads, msg));
        double ttot = 0, ttrial[100];
        for (int trial = 0; trial < ntrials; trial++) {
            tt = LAGraph_WallClockTime();
            LAGRAPH_TRY (LAGraph_cdlp(&c, G, 100, msg));
            GRB_TRY (GrB_free(&c));
            ttrial[trial] = LAGraph_WallClockTime() - tt;
            ttot += ttrial[trial];
            printf("threads %2d trial %2d: %12.6f sec\n",
                   nthreads, trial, ttrial[trial]);
            fprintf(stderr, "threads %2d trial %2d: %12.6f sec\n",
                    nthreads, trial, ttrial[trial]);
        }
        ttot = ttot / ntrials;

        printf("Avg: CDLP "
               "nthreads: %3d time: %12.6f matrix: %s\n",
               nthreads, ttot, matrix_name);

        fprintf(stderr, "Avg: CDLP "
                        "nthreads: %3d time: %12.6f matrix: %s\n",
                nthreads, ttot, matrix_name);
    }

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
