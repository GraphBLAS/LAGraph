//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/matching_demo.c: benchmarks for
// LAGr_MaximumMatching
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

// #define VERBOSE

#define NTHREAD_LIST 1
#define THREAD_LIST 0

#define LG_FREE_ALL                                                            \
    {                                                                          \
        GrB_free(&A);                                                          \
        GrB_free(&M);                                                          \
        GrB_free(&mateC);                                                      \
        LAGraph_Free((void *)&I, NULL);                                        \
        LAGraph_Free((void *)&X, NULL);                                        \
    }

#if LAGRAPH_SUITESPARSE

GrB_Info check_matching(GrB_Matrix A, GrB_Vector mateC, char *msg)
{
    GrB_Index nmatched = 0;
    GrB_Vector mateR = NULL;
    GrB_Matrix M = NULL;
    GrB_Index *I = NULL, *X = NULL;

    uint64_t ncols = 0, nrows = 0;
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));

    // invert to check for dups
    GrB_Index IBytes = 0, XBytes = 0;
    bool jumbled = 1;
    GRB_TRY(GrB_Vector_new(&mateR, GrB_UINT64, nrows));
    GRB_TRY(GxB_Vector_unpack_CSC(mateC, (GrB_Index **)&I, (void **)&X, &IBytes,
                                  &XBytes, NULL, &nmatched, &jumbled, NULL));
    GRB_TRY(GrB_Vector_build_UINT64(mateR, X, I, nmatched, GrB_FIRST_UINT64));
    GrB_Index nmateR = 0;
    GRB_TRY(GrB_Vector_nvals(&nmateR, mateR));
    // if nvals of mateC and mateR don't match, then there's at least
    // one row that is used in at least one matching
    if (nmatched != nmateR)
    {
        printf("Duplicates in mateC");
        fflush(stdout);
        abort();
    }

    // pack matched values in a matrix
    bool *val;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&val, nmatched, sizeof(bool), msg));
    for (uint64_t i = 0; i < nmatched; i++)
        val[i] = 1;
    GRB_TRY(GrB_Matrix_new(&M, GrB_BOOL, nrows, ncols));
    GRB_TRY(GrB_Matrix_build_BOOL(M, X, I, val, nmatched, NULL));
    LAGRAPH_TRY(LAGraph_Free((void **)&val, msg));
    // mask with matrix A to check if all edges are present in A
    GRB_TRY(GrB_Matrix_assign(M, M, NULL, A, GrB_ALL, nrows, GrB_ALL, ncols,
                              GrB_DESC_S));
    GrB_Index nvalsM = 0;
    GRB_TRY(GrB_Matrix_nvals(&nvalsM, M));
    // if values have been eliminated then edges do not exist in A
    if (nvalsM != nmatched)
    {
        printf("mateC invalid!\n");
        fflush(stdout);
        abort();
    }

    GRB_TRY(GxB_Vector_pack_CSC(mateC, (GrB_Index **)&I, (void **)&X, IBytes,
                                XBytes, false, nmatched, jumbled, NULL));

    GrB_Vector_free(&mateR);
    GrB_Matrix_free(&M);
    return (GrB_SUCCESS);
}
#endif

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
#if LAGRAPH_SUITESPARSE

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

    if (argc < 2)
    {
        printf("Invalid usage, please read comments\n");
        fflush(stdout);
        return 0;
    }
    char *matrix_name = (argc > 1) ? argv[1] : "stdin";

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
    fflush(stdout);
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
    fflush(stdout);
#endif

    //--------------------------------------------------------------------------
    // warmup before benchmarking
    //--------------------------------------------------------------------------

    double t = LAGraph_WallClockTime();
    LAGRAPH_TRY(
        LAGr_MaximumMatching(&mateC, NULL, A, AT, mateC_init, true, msg));
    t = LAGraph_WallClockTime() - t;
    LAGRAPH_TRY(check_matching(A, mateC, msg));
    uint64_t sprank = 0;
    GRB_TRY(GrB_Vector_nvals(&sprank, mateC));
    printf("number of matches: %" PRIu64 "\n", sprank);
    fflush(stdout);
    GRB_TRY(GrB_free(&mateC));
#ifdef VERBOSE
    printf("warmup time %g sec\n", t);
    fflush(stdout);
#endif

    //--------------------------------------------------------------------------
    // benchmark
    //--------------------------------------------------------------------------

    // the GAP benchmark requires 16 trials
    int ntrials = 16;
#ifdef VERBOSE
    printf("# of trials: %d\n", ntrials);
    fflush(stdout);
#endif

    for (int kk = 1; kk <= nt; kk++)
    {
        int nthreads = Nthreads[kk];
        if (nthreads > nthreads_max)
            continue;
        LAGRAPH_TRY(LAGraph_SetNumThreads(1, nthreads, msg));

#ifdef VERBOSE
        printf("\n--------------------------- nthreads: %2d\n", nthreads);
        fflush(stdout);
#endif

        double total_time = 0;

        for (int trial = 0; trial < ntrials; trial++)
        {
            t = LAGraph_WallClockTime();
            LAGRAPH_TRY(LAGr_MaximumMatching(&mateC, NULL, A, AT, mateC_init,
                                             true, msg));
            t = LAGraph_WallClockTime() - t;
            GRB_TRY(GrB_free(&mateC));
#ifdef VERBOSE
            printf("trial: %2d time: %10.7f sec\n", trial, t);
            fflush(stdout);
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
        fflush(stdout);
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------
    LG_FREE_ALL;

    LAGRAPH_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
