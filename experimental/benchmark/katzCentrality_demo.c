//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/katzCentrality_demo.c: benchmark for Katz
//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include <stdlib.h>

#define NTHREAD_LIST 1
#define THREAD_LIST 0

#define LG_FREE_ALL             \
{                               \
    LAGraph_Delete (&G, NULL) ; \
    GrB_free (&c) ;             \
}

int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    LAGraph_Graph G = NULL ;
    GrB_Vector c = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    int ntrials = 3 ;
    printf ("# of trials: %d\n", ntrials) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;

    int nthreads_max, nthreads_outer, nthreads_inner ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
    nthreads_max = nthreads_outer * nthreads_inner ;

    if (Nthreads [1] == 0)
    {
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0)
                nt = t-1 ;
        }
    }

    printf ("threads to test:") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max)
            continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // read in the graph and parse optional Katz parameters
    //--------------------------------------------------------------------------

    // Usage:
    //   katzCentrality_demo <matrix-market-file> <alpha>
    //                       [beta] [max_iters] [tol] [normalize] [use_weights]
    // where normalize and use_weights are 0/1.
    if (argc != 1 && (argc < 3 || argc > 8))
    {
        printf ("Usage: %s <matrix-market-file> <alpha>"
                " [beta] [max_iters] [tol] [normalize] [use_weights]\n",
                argv [0]) ;
        return (GrB_INVALID_VALUE) ;
    }

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;

    double  alpha      = 0.01 ;
    double  beta       = 1.00 ;
    int64_t iters      = 0 ;
    int64_t max_iter   = 1000 ;
    double  tol        = 1e-6 ;
    bool    normalize  = false ;
    bool    use_weights = false ;

    if (argc > 2) alpha       = strtod (argv [2], NULL) ;
    if (argc > 3) beta        = strtod (argv [3], NULL) ;
    if (argc > 4) max_iter    = (int64_t) strtoll (argv [4], NULL, 10) ;
    if (argc > 5) tol         = strtod (argv [5], NULL) ;
    if (argc > 6) normalize   = (atoi (argv [6]) != 0) ;
    if (argc > 7) use_weights = (atoi (argv [7]) != 0) ;

    LAGRAPH_TRY(
        readproblem(&G, NULL, false, false, false, NULL, false, argc, argv));

    // Katz uses incoming edges; AT is required for directed graphs.
    int cache_result = LAGraph_Cached_AT (G, msg) ;
    LG_ASSERT_MSG (cache_result == GrB_SUCCESS ||
                   cache_result == LAGRAPH_CACHE_NOT_NEEDED,
                   cache_result, "LAGraph_Cached_AT failed") ;

    //--------------------------------------------------------------------------
    // benchmark Katz centrality
    //--------------------------------------------------------------------------

    printf ("\nKatz params: alpha=%g beta=%g max_iter=%lld tol=%g"
            " normalize=%d use_weights=%d\n\n",
            alpha, beta, (long long) max_iter, tol,
            (int) normalize, (int) use_weights) ;

    // warmup for more accurate timing
    double tt = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGr_KatzCentrality (&c, &iters, G, alpha, beta,
        max_iter, tol, normalize, use_weights, msg)) ;
    tt = LAGraph_WallClockTime () - tt ;
    GRB_TRY (GrB_free (&c)) ;
    printf ("warmup time %g sec\n", tt) ;

    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max)
            continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;

        double ttot = 0, ttrial [100] ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            double t1 = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGr_KatzCentrality (&c, &iters, G, alpha, beta,
                max_iter, tol, normalize, use_weights, msg)) ;
            GRB_TRY (GrB_free (&c)) ;
            ttrial [trial] = LAGraph_WallClockTime () - t1 ;
            ttot += ttrial [trial] ;

            printf ("threads %2d trial %2d: %12.6f sec\n",
                    nthreads, trial, ttrial [trial]) ;
            fprintf (stderr, "threads %2d trial %2d: %12.6f sec\n",
                     nthreads, trial, ttrial [trial]) ;
        }

        ttot = ttot / ntrials ;
        printf (         "Avg: KatzCentrality nthreads: %3d time: %12.6f sec"
                         " (%" PRId64 " iterations) matrix: %s\n",
                         nthreads, ttot, iters, matrix_name) ;
        fprintf (stderr, "Avg: KatzCentrality nthreads: %3d time: %12.6f sec"
                         " (%" PRId64 " iterations) matrix: %s\n",
                         nthreads, ttot, iters, matrix_name) ;
    }

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
