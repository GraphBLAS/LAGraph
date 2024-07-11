//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/matching_demo.c: benchmarks for
// LAGraph_MaximumMatching
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Christina Koutsou, Aristotle University of Thessaloniki

/*
Usage:

*/

//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"
#include <omp.h>

// #define VERBOSE

#define NTHREAD_LIST 1
#define THREAD_LIST 0

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                            \
    {                                                                          \
        LAGraph_Delete(&G, NULL);                                              \
        GrB_free(&A);                                                          \
        GrB_free(&AT);                                                         \
        GrB_free(&mateC_init);                                                 \
        GrB_free(&mateC);                                                      \
    }

int main(int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // declare inputs and outputs
    //--------------------------------------------------------------------------

    char msg[LAGRAPH_MSG_LEN];

    LAGraph_Graph G = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix AT = NULL;
    GrB_Vector mateC_init = NULL;
    GrB_Vector mateC = NULL;

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    bool burble = false;
    demo_init(burble);

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    // if (argc < 2)
    // {
    //     printf("Invalid usage, please read comments\n");
    //     return 0;
    // }
    char *matrix_name = (argc > 1) ? argv[1] : "stdin";

    LAGRAPH_TRY(LAGraph_Random_Init(msg));
    bool make_symmetric = false, remove_self_edges = false, structural = true,
         ensure_positive = false;
    LAGRAPH_TRY(readproblem(&G, NULL, make_symmetric, remove_self_edges,
                            structural, NULL, ensure_positive, argc, argv));

    A = G->A;
    // compute AT to be able to use push-pull optimization
    if (G->is_symmetric_structure)
        AT = A;
    else
    {
        LAGRAPH_TRY(LAGraph_Cached_AT(G, msg));
        AT = G->AT;
    }

    //--------------------------------------------------------------------------
    // determine the number of threads to run the algorithm with
    //--------------------------------------------------------------------------

    int nt = NTHREAD_LIST;
    int Nthreads[20] = {0, THREAD_LIST};
    int nthreads_max, nthreads_outer, nthreads_inner;
    LAGRAPH_TRY(LAGraph_GetNumThreads(&nthreads_outer, &nthreads_inner, msg));
#ifdef VERBOSE
    printf("nthreads_outer: %d, nthreads_inner: %d\n", nthreads_outer,
           nthreads_inner);
#endif
    nthreads_max = nthreads_outer * nthreads_inner;
    if (Nthreads[1] == 0) // THREAD_LIST == 0
    {
        // create thread list automatically
        Nthreads[1] = nthreads_max;
        for (int t = 2; t <= nt; t++)
        {
            Nthreads[t] = Nthreads[t - 1] / 2;
            if (Nthreads[t] == 0)
                nt = t - 1;
        }
    }
#ifdef VERBOSE
    printf("threads to test: ");
    for (int t = 1; t <= nt; t++)
    {
        int nthreads = Nthreads[t];
        if (nthreads > nthreads_max)
            continue;
        printf(" %d", nthreads);
    }
    printf("\n");
#endif

    //--------------------------------------------------------------------------
    // warmup before benchmarking
    //--------------------------------------------------------------------------

    double t = LAGraph_WallClockTime();
    LAGRAPH_TRY(LAGraph_MaximumMatching(&mateC, A, AT, mateC_init, msg));
    t = LAGraph_WallClockTime() - t;
    GRB_TRY(GrB_free(&mateC));
#ifdef VERBOSE
    printf("warmup time %g sec\n", t);
#endif

    //--------------------------------------------------------------------------
    // benchmark
    //--------------------------------------------------------------------------

    // the GAP benchmark requires 16 trials
    int ntrials = 16;
    // ntrials = 1 ;    // HACK to run just one trial
#ifdef VERBOSE
    printf("# of trials: %d\n", ntrials);
#endif

    for (int kk = 1; kk <= nt; kk++)
    {
        int nthreads = Nthreads[kk];
        if (nthreads > nthreads_max)
            continue;
        LAGRAPH_TRY(LAGraph_SetNumThreads(1, nthreads, msg));

#ifdef VERBOSE
        printf("\n--------------------------- nthreads: %2d\n", nthreads);
#endif

        double total_time = 0;

        for (int trial = 0; trial < ntrials; trial++)
        {
            t = LAGraph_WallClockTime();
            LAGRAPH_TRY(
                LAGraph_MaximumMatching(&mateC, A, AT, mateC_init, msg));
            t = LAGraph_WallClockTime() - t;
            GRB_TRY(GrB_free(&mateC));
#ifdef VERBOSE
            printf("trial: %2d time: %10.7f sec\n", trial, t);
#endif
            total_time += t;
        }

        double total_time_per_trial = total_time / ntrials;

#ifndef VERBOSE
        printf("%.7f\n", total_time_per_trial);
#else
        printf("maximum matching: %3d: avg time: %10.7f (sec) matrix: %s\n",
               nthreads, total_time_per_trial, matrix_name);
#endif
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------
    LG_FREE_ALL;

    LAGRAPH_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS);
}
