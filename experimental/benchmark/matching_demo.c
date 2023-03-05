#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"
#include <omp.h>

// #define VERBOSE

#define NTHREAD_LIST 4
#define THREAD_LIST 40, 20, 16, 8, // 4, 2, 1

int main (int argc, char** argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    bool test_performance = true ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix E = NULL ;
    GrB_Vector matching = NULL ;

    bool burble = false ; 
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    int force_stdin = 0 ;

    if (argc > 1) {
        // stdin is never a valid filename; if this is in argv[1], then it means I want to use this demo run
        // to print out data for my external tests, and not benchmark performance; I have my own arguments stored in argv.
        force_stdin = ( strcmp (argv [1], "stdin") == 0 ) ;

        if (force_stdin) {
            // mark that I am not running performance benchmarks, but printing data for my external tests
            test_performance = false ;
        }
    }
    LAGRAPH_TRY (LAGraph_Random_Init (msg)) ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, false, GrB_FP64, false, force_stdin ? 1 : argc, argv)) ;

    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    GRB_TRY (LAGraph_A_to_E (&E, G, msg)) ;

    if (!test_performance) {

        //--------------------------------------------------------------------------
        // Printing E matrix, best result from ntrial runs for my own, external tests for quality (not performance)
        //--------------------------------------------------------------------------
        int ntrials = atoi(argv [3]) ;

        GrB_Vector best_matching = NULL ;
        GrB_Index max_val = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++) {
            int64_t seed = trial * n + 1 ;
            LAGRAPH_TRY (LAGraph_MaximalMatching (&matching, E, atoi(argv [2]), seed, msg)) ;
            GrB_Index nvals ;
            GRB_TRY (GrB_Vector_nvals (&nvals, matching)) ;
            if (nvals > max_val) {
                if (best_matching != NULL) {
                    GrB_free (&best_matching) ;
                }
                best_matching = matching ;
                max_val = nvals ;
            } else {
                GrB_free (&matching) ;
            }
        }
        LAGRAPH_TRY (LAGraph_Vector_Print (best_matching, LAGraph_COMPLETE, stdout, msg)) ;
        LAGRAPH_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;

        GrB_free (&best_matching) ;
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

        printf("%.7f\n", t);

#ifdef VERBOSE
        printf ("weighted max matching: %3d: avg time: %10.7f (sec) matrix: %s\n",
                nthreads, t, matrix_name) ;
        fprintf (stderr, "weighted max matching: %3d: avg time: %10.7f (sec) matrix: %s\n",
                nthreads, t, matrix_name) ;
#endif
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------
    LAGraph_Delete (&G, NULL) ;
    GrB_free (&E) ;

    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}