//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/katzCentrality_demo.c: benchmark for Katz
//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0

#define LG_FREE_ALL               \
    {                             \
        LAGraph_Delete(&G, NULL); \
        GrB_free(&c);             \
    }

int main(int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg[LAGRAPH_MSG_LEN];
    LAGraph_Graph G = NULL;
    GrB_Vector c = NULL;

    // start GraphBLAS and LAGraph
    bool burble = false;
    demo_init(burble);

    int ntrials = 3;
    printf("# of trials: %d\n", ntrials);

    int nt = NTHREAD_LIST;
    int Nthreads[20] = {0, THREAD_LIST};

    int nthreads_max, nthreads_outer, nthreads_inner;
    LAGRAPH_TRY(LAGraph_GetNumThreads(&nthreads_outer, &nthreads_inner, msg));
    nthreads_max = nthreads_outer * nthreads_inner;

    if (Nthreads[1] == 0)
    {
        Nthreads[1] = nthreads_max;
        for (int t = 2; t <= nt; t++)
        {
            Nthreads[t] = Nthreads[t - 1] / 2;
            if (Nthreads[t] == 0)
                nt = t - 1;
        }
    }

    printf("threads to test:");
    for (int t = 1; t <= nt; t++)
    {
        int nthreads = Nthreads[t];
        if (nthreads > nthreads_max)
            continue;
        printf(" %d", nthreads);
    }
    printf("\n");

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv[1] : "stdin";
    LAGRAPH_TRY(readproblem(&G, NULL,
                            true, true, true, NULL, false, argc, argv));

    // Katz uses incoming edges; AT is required for directed graphs.
    int cache_result = LAGraph_Cached_AT(G, msg);
    LG_ASSERT_MSG(cache_result == GrB_SUCCESS ||
                      cache_result == LAGRAPH_CACHE_NOT_NEEDED,
                  cache_result, "LAGraph_Cached_AT failed");

    //--------------------------------------------------------------------------
    // benchmark Katz centrality
    //--------------------------------------------------------------------------

    double alpha = 0.01;
    double beta = 1.00;
    int max_iter = 1000;
    double tol = 1e-6;
    bool normalize = false;

    // warmup for more accurate timing
    double tt = LAGraph_WallClockTime();
    LAGRAPH_TRY(LAGr_KatzCentrality(&c, G, alpha, beta,
                                    max_iter, tol, normalize, msg));
    tt = LAGraph_WallClockTime() - tt;
    GRB_TRY(GrB_free(&c));
    printf("warmup time %g sec\n", tt);

    for (int t = 1; t <= nt; t++)
    {
        int nthreads = Nthreads[t];
        if (nthreads > nthreads_max)
            continue;
        LAGRAPH_TRY(LAGraph_SetNumThreads(1, nthreads, msg));

        double ttot = 0, ttrial[100];
        for (int trial = 0; trial < ntrials; trial++)
        {
            double t1 = LAGraph_WallClockTime();
            LAGRAPH_TRY(LAGr_KatzCentrality(&c, G, alpha, beta,
                                            max_iter, tol, normalize, msg));
            GRB_TRY(GrB_free(&c));
            ttrial[trial] = LAGraph_WallClockTime() - t1;
            ttot += ttrial[trial];

            printf("threads %2d trial %2d: %12.6f sec\n",
                   nthreads, trial, ttrial[trial]);
            fprintf(stderr, "threads %2d trial %2d: %12.6f sec\n",
                    nthreads, trial, ttrial[trial]);
        }

        ttot = ttot / ntrials;
        printf("Avg: KatzCentrality nthreads: %3d time: %12.6f matrix: %s\n",
               nthreads, ttot, matrix_name);
        fprintf(stderr, "Avg: KatzCentrality nthreads: %3d time: %12.6f matrix: %s\n",
                nthreads, ttot, matrix_name);
    }

    LG_FREE_ALL;
    LAGRAPH_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS);
}
