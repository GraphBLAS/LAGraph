//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/matching_demo.c: runner for LAGraph_MaximalMatching
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Vidith Madhu, Texas A&M University

/*
Usage:
Option 1: Run for performance
./matching_demo <matrix_name> <matching_type>
matrix_name: either the name of the .mtx file or empty for stdin
matching_type: 0, 1, 2 for random matching, heavy edge matching, and light edge matching respectively
NOTE: This is the typical scenario that all other benchmark codes are used for

Option 2: Run for quality
./matching_demo -q <matching_type> <ntrials>
NOTE: this option only accepts input via stdin
-q option as the matrix name specifies to run for quality, not performance
matching_type: 0, 1, 2 for random matching, heavy edge matching, and light edge matching respectively
ntrials: How many trials to run (picks matching with the highest quality across all runs)
NOTE: When complete, prints out the matching vector and E matrix of the input graph. This can be piped out to any file.
*/

//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"
#include <omp.h>

#define VERBOSE

#define NTHREAD_LIST 1
#define THREAD_LIST 8

#undef LG_FREE_ALL
#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&E) ;                 \
    GrB_free (&E_t) ;               \
    GrB_free (&matching) ;          \
    GrB_free (&weight) ;            \
    GrB_free (&best_matching) ;     \
    GrB_free (&use_weights) ;       \
}

int main (int argc, char** argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    bool test_performance = true ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix E = NULL ;
    GrB_Matrix E_t = NULL ;
    GrB_Vector matching = NULL ;
    GrB_Vector weight = NULL ;
    GrB_Vector best_matching = NULL ;
    GrB_Vector use_weights = NULL ;

    bool burble = false ; 
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    int force_stdin = 0 ;

    if (argc > 1) {
        // -q option as the matrix name means to run the quality tests
        force_stdin = ( strcmp (argv [1], "-q") == 0 ) ;

        if (force_stdin) {
            // mark that I am not running performance benchmarks, but printing data for my external tests
            test_performance = false ;
        }
    }
    LAGRAPH_TRY (LAGraph_Random_Init (msg)) ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, false, GrB_FP64, false, force_stdin ? 1 : argc, argv)) ;

    GrB_Index n ;
    GrB_Index num_edges ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    GRB_TRY (LAGraph_Incidence_Matrix (&E, G, msg)) ;
    GRB_TRY (GrB_Matrix_ncols (&num_edges, E)) ;

    GRB_TRY (GrB_Matrix_new (&E_t, GrB_FP64, num_edges, n)) ;
    GRB_TRY (GrB_Vector_new (&weight, GrB_FP64, num_edges)) ;
    
    GRB_TRY (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;

    GRB_TRY (GrB_reduce (weight, NULL, NULL, GrB_MAX_MONOID_FP64, E_t, NULL)) ;

    if (!test_performance) {

        //--------------------------------------------------------------------------
        // Printing E matrix, best result from ntrial runs for my own, external tests for quality (not performance)
        //--------------------------------------------------------------------------
        int ntrials = atoi(argv [3]) ;
        int matching_type = atoi(argv [2]) ;

        // best answer so far
        double best_val = (matching_type == 2) ? 1e18 : 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++) {
            int64_t seed = trial * n + 1 ;
            
            LAGRAPH_TRY (LAGraph_MaximalMatching (&matching, E, matching_type, seed, msg)) ;
            double matching_value = 0 ;
            if (matching_type != 0) {
                // weighted matching; need to compute total weight of matching
                GRB_TRY (GrB_Vector_new (&use_weights, GrB_FP64, num_edges)) ;
                GRB_TRY (GrB_eWiseMult (use_weights, NULL, NULL, GrB_TIMES_FP64, weight, matching, NULL)) ;
                GRB_TRY (GrB_reduce (&matching_value, NULL, GrB_PLUS_MONOID_FP64, use_weights, NULL)) ;
                GRB_TRY (GrB_free (&use_weights)) ;
            } else {
                // random matching; just want number of matched edges
                GrB_Index matching_value_int ;
                GRB_TRY (GrB_Vector_nvals (&matching_value_int, matching)) ;
                matching_value = matching_value_int ;
            }
            bool cond = (matching_value > best_val) ;
            if (matching_type == 2) {
                // for light matchings, want to prefer smaller values
                cond = (matching_value < best_val) ;
            }
            if (cond) {
                if (best_matching != NULL) {
                    // free the previous vector before updating
                    GRB_TRY (GrB_free (&best_matching)) ;
                }
                best_matching = matching ;
                best_val = matching_value ;
            } else {
                GRB_TRY (GrB_free (&matching)) ;
            }
        }
        // print matching vector and E matrix; matching will be validated again in custom tests
        LAGRAPH_TRY (LAGraph_Vector_Print (best_matching, LAGraph_COMPLETE, stdout, msg)) ;
        LAGRAPH_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;

        LG_FREE_ALL ;
        
        return 0 ;
    }
    int nt = NTHREAD_LIST ;
    
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max, nthreads_outer, nthreads_inner ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, NULL)) ;
#ifdef VERBOSE
    printf("nthreads_outer: %d, nthreads_inner: %d\n", nthreads_outer, nthreads_inner);
#endif
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
#ifdef VERBOSE
    printf ("threads to test: ") ;
#endif
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
#ifdef VERBOSE
        printf (" %d", nthreads) ;
#endif
    }
#ifdef VERBOSE
    printf ("\n") ;
#endif

    // warmup for more accurate timing
    double tt = LAGraph_WallClockTime ( ) ;
    // user-provided matching type (random, heavy, light)
    int match_type = (argc > 2) ? atoi(argv[2]) : 1;
    // GRB_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;
    LAGRAPH_TRY (LAGraph_MaximalMatching (&matching, E, match_type, 5, msg)) ;
    tt = LAGraph_WallClockTime ( ) - tt ;
    GRB_TRY (GrB_free (&matching)) ;
#ifdef VERBOSE
    printf ("warmup time %g sec\n", tt) ;
#endif
    GrB_free (&matching) ;

    // the GAP benchmark requires 16 trials
    int ntrials = 16 ;
    // ntrials = 1 ;    // HACK to run just one trial
#ifdef VERBOSE
    printf ("# of trials: %d\n", ntrials) ;
#endif

    for (int kk = 1 ; kk <= nt ; kk++)
    {
        int nthreads = Nthreads [kk] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;

#ifdef VERBOSE
        printf ("\n--------------------------- nthreads: %2d\n", nthreads) ;
#endif

        double total_time = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t seed = trial * n + 1 ;
            double tt = LAGraph_WallClockTime ( ) ;
            LAGRAPH_TRY (LAGraph_MaximalMatching (&matching, E, match_type, seed, msg)) ;
            tt = LAGraph_WallClockTime ( ) - tt ;
            GRB_TRY (GrB_free (&matching)) ;
#ifdef VERBOSE
            printf ("trial: %2d time: %10.7f sec\n", trial, tt) ;
#endif
            total_time += tt ;
        }

        double t = total_time / ntrials ;

#ifndef VERBOSE
        printf("%.7f\n", t);
#endif

#ifdef VERBOSE
        printf ("maximal matching: %3d: avg time: %10.7f (sec) matrix: %s\n",
                nthreads, t, matrix_name) ;
        fprintf (stderr, "maximal matching: %3d: avg time: %10.7f (sec) matrix: %s\n",
                nthreads, t, matrix_name) ;
#endif
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------
    LG_FREE_ALL ;

    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
