//------------------------------------------------------------------------------
// walks_demo.c: benchmark LAGraph_NumberOfWalks vs BFS algorithms
//------------------------------------------------------------------------------

//  Benchmarks two parts of the NumberOfWalks algorithm:
//
//  All-pairs:
//      LAGraph_NumberOfWalks vs LAGraph_MultiSourceBFS
//      Both produce an n x (n or nsrc) matrix.
//      NumberOfWalks uses binary exponentiation (O(log k) mxm calls)
//      MultiSourceBFS uses k mxm calls.
//
//  Single-source:
//      LAGraph_NumberOfWalks (src=indicator) vs LAGr_BreadthFirstSearch
//      Both start from one node and produce a length-n result vector/row.
//
//  Usage:
//      ./walks_demo graph.mtx                      (generates default sources)
//      ./walks_demo graph.mtx sources.mtx          (sources in GAP format)
//      ./walks_demo graph.mtx sources.mtx 6        (walk length k=6)
//
// The sources file is a Matrix Market file with one source node index
// per row in column 0, matching the GAP benchmark format used by bfs_demo. 
// If omitted, the first N_DEFAULT_SOURCES node indices are used.

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0
#define DEFAULT_K         4
#define N_DEFAULT_SOURCES 64
#define MAX_TRIALS        8

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                     \
{                                                       \
    LAGraph_Delete (&G, msg) ;                          \
    GrB_free (&C) ;                                     \
    GrB_free (&level) ;                                 \
    GrB_free (&parent) ;                                \
    GrB_free (&src_indicator) ;                         \
    GrB_free (&multisrc_vector) ;                       \
    GrB_free (&SourceNodes) ;                           \
    LAGraph_Free ((void **) &t_walks_all,  msg) ;       \
    LAGraph_Free ((void **) &t_msbfs,      msg) ;       \
    LAGraph_Free ((void **) &t_walks_src,  msg) ;       \
    LAGraph_Free ((void **) &t_bfs,        msg) ;       \
}

int main (int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix C = NULL ;
    GrB_Matrix level = NULL, parent = NULL ;
    GrB_Vector src_indicator  = NULL ;  // indicator for SS walks
    GrB_Vector multisrc_vector = NULL ; // source list for MultiSourceBFS
    GrB_Matrix SourceNodes = NULL ;

    double *t_walks_all = NULL ;        // NumberOfWalks time per thread count
    double *t_msbfs     = NULL ;        // MultiSourceBFS time per thread count

    double *t_walks_src = NULL ;        // avg SS NumberOfWalks per thread count
    double *t_bfs       = NULL ;        // avg BFS per thread count

    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read walk length k from the WALK_K environment variable
    //--------------------------------------------------------------------------

    int64_t k = DEFAULT_K ;
    const char *k_env = getenv ("WALK_K") ;
    if (k_env != NULL && k_env [0] != '\0') k = atol (k_env) ;
    if (k <= 0) k = DEFAULT_K ;

    //--------------------------------------------------------------------------
    // THREAD SETUP
    //--------------------------------------------------------------------------

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
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test:") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_walks_all, nthreads_max+1,
        sizeof (double), msg)) ;
    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_msbfs,     nthreads_max+1,
        sizeof (double), msg)) ;
    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_walks_src, nthreads_max+1,
        sizeof (double), msg)) ;
    LAGRAPH_TRY (LAGraph_Malloc ((void **) &t_bfs,       nthreads_max+1,
        sizeof (double), msg)) ;

    //--------------------------------------------------------------------------
    // LOAD
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, &SourceNodes,
        false,      // make_symmetric: walks work on directed graphs
        false,      // remove_self_edges: self loops affect walk counts
        false,      // structural: need INT64 values, not bool
        GrB_INT64,  // typecast all entries to INT64
        false,      // ensure_positive: walks can be zero
        argc, argv)) ;

    LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // build source node list
    // If no source file was given, use the first N_DEFAULT_SOURCES nodes.
    //--------------------------------------------------------------------------

    GrB_Index ntrials ;
    if (SourceNodes == NULL)
    {
        ntrials = (GrB_Index) N_DEFAULT_SOURCES ;
        if (ntrials > n) ntrials = n ;
        GRB_TRY (GrB_Matrix_new (&SourceNodes, GrB_INT64, ntrials, 1)) ;
        for (GrB_Index i = 0 ; i < ntrials ; i++)
        {
            GRB_TRY (GrB_Matrix_setElement_INT64 (SourceNodes,
                (int64_t)(i + 1), i, 0)) ;
        }
    }
    else
    {
        GRB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;
    }
    if (ntrials > MAX_TRIALS) ntrials = MAX_TRIALS ;

    printf ("graph: %s  n: %" PRIu64 "  edges: %" PRIu64 "  k: %" PRId64 "\n",
        matrix_name, (uint64_t) n, (uint64_t) nvals, k) ;
    printf ("source trials: %" PRIu64 "\n", (uint64_t) ntrials) ;
    fflush (stdout) ; fflush (stderr) ;

    //--------------------------------------------------------------------------
    // build MultiSourceBFS source vector
    // GrB_Vector where entry i = source node index (0-based)
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&multisrc_vector, GrB_INT64, ntrials)) ;
    for (GrB_Index i = 0 ; i < ntrials ; i++)
    {
        int64_t src ;
        GRB_TRY (GrB_Matrix_extractElement_INT64 (&src, SourceNodes, i, 0)) ;
        src-- ;  // convert 1-based (GAP format) to 0-based
        GRB_TRY (GrB_Vector_setElement_INT64 (multisrc_vector, src, i)) ;
    }

    //--------------------------------------------------------------------------
    // WARMUP
    //--------------------------------------------------------------------------

    printf ("\n--- warmup ---\n") ;

    // warmup: all pairs NumberOfWalks
    double twarm = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, NULL, k)) ;
    twarm = LAGraph_WallClockTime () - twarm ;
    GRB_TRY (GrB_free (&C)) ;
    printf ("NumberOfWalks (all-pairs): %g sec\n", twarm) ;

    // warmup: MultiSourceBFS
    twarm = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, NULL, G, multisrc_vector, msg)) ;
    twarm = LAGraph_WallClockTime () - twarm ;
    GRB_TRY (GrB_free (&level)) ;
    printf ("MultiSourceBFS:            %g sec\n", twarm) ;

    // warmup: SS BFS and NumberOfWalks
    {
        int64_t src ;
        GRB_TRY (GrB_Matrix_extractElement_INT64 (&src, SourceNodes, 0, 0)) ;
        src-- ;

        twarm = LAGraph_WallClockTime () ;
        LAGRAPH_TRY (LAGr_BreadthFirstSearch (NULL, (GrB_Vector *) &parent, G,
            (GrB_Index) src, msg)) ;
        twarm = LAGraph_WallClockTime () - twarm ;
        GRB_TRY (GrB_free (&parent)) ;
        printf ("BFS (single-source):       %g sec\n", twarm) ;

        GRB_TRY (GrB_Vector_new (&src_indicator, GrB_INT64, n)) ;
        GRB_TRY (GrB_Vector_setElement_INT64 (src_indicator, 1, (GrB_Index) src)) ;
        twarm = LAGraph_WallClockTime () ;
        LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, src_indicator, k)) ;
        twarm = LAGraph_WallClockTime () - twarm ;
        GRB_TRY (GrB_free (&C)) ;
        GRB_TRY (GrB_free (&src_indicator)) ;
        printf ("NumberOfWalks (single-source): %g sec\n", twarm) ;
    }
    fflush (stdout) ; fflush (stderr) ;

    //==========================================================================
    // BENCHMARK
    //==========================================================================

    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;

        printf ("\n=========================== nthreads: %2d ===========================\n",
            nthreads) ;

        //----------------------------------------------------------------------
        // All pairs NumberOfWalks vs MultiSourceBFS
        //----------------------------------------------------------------------

        printf ("\n--- Part 1: All-pairs (k=%" PRId64 ", %" PRIu64 " sources) ---\n",
            k, (uint64_t) ntrials) ;

        // NumberOfWalks: produces the full A^k matrix
        {
            double t_run = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, NULL, k)) ;
            t_run = LAGraph_WallClockTime () - t_run ;
            GRB_TRY (GrB_free (&C)) ;
            t_walks_all [nthreads] = t_run ;
            printf ("NumberOfWalks all-pairs  k: %2" PRId64
                "  threads: %2d  time: %10.4f sec\n",
                k, nthreads, t_run) ;
            fflush (stdout) ;
        }

        // MultiSourceBFS: over all ntrials source nodes
        {
            double t_run = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, NULL, G,
                multisrc_vector, msg)) ;
            t_run = LAGraph_WallClockTime () - t_run ;
            GRB_TRY (GrB_free (&level)) ;
            t_msbfs [nthreads] = t_run ;
            printf ("MultiSourceBFS           sources: %" PRIu64
                "  threads: %2d  time: %10.4f sec\n",
                (uint64_t) ntrials, nthreads, t_run) ;
            fflush (stdout) ;
        }

        //----------------------------------------------------------------------
        // SS NumberOfWalks vs BFS, averaged over all sources
        //----------------------------------------------------------------------

        printf ("\n--- Part 2: Single-source (k=%" PRId64 ") ---\n", k) ;

        double total_walks_src = 0, total_bfs = 0 ;
        GRB_TRY (GrB_Vector_new (&src_indicator, GrB_INT64, n)) ;

        for (GrB_Index trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t src ;
            GRB_TRY (GrB_Matrix_extractElement_INT64 (&src, SourceNodes,
                trial, 0)) ;
            src-- ;

            // SS BFS
            double tb = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGr_BreadthFirstSearch (NULL, (GrB_Vector *) &parent,
                G, (GrB_Index) src, msg)) ;
            tb = LAGraph_WallClockTime () - tb ;
            GRB_TRY (GrB_free (&parent)) ;
            total_bfs += tb ;

            // SS NumberOfWalks
            GRB_TRY (GrB_Vector_setElement_INT64 (src_indicator, 1,
                (GrB_Index) src)) ;
            double tw = LAGraph_WallClockTime () ;
            LAGRAPH_TRY (LAGraph_NumberOfWalks (&C, G->A, src_indicator, k)) ;
            tw = LAGraph_WallClockTime () - tw ;
            GRB_TRY (GrB_free (&C)) ;
            // clear for next trial
            GRB_TRY (GrB_Vector_clear (src_indicator)) ;
            total_walks_src += tw ;

            printf ("trial: %3" PRIu64 "  src: %12" PRId64
                "  walks: %8.4f sec  bfs: %8.4f sec\n",
                (uint64_t) trial, src, tw, tb) ;
            fflush (stdout) ;
        }

        GRB_TRY (GrB_free (&src_indicator)) ;

        t_walks_src [nthreads] = total_walks_src / (double) ntrials ;
        t_bfs       [nthreads] = total_bfs       / (double) ntrials ;

        //----------------------------------------------------------------------
        // SUMMARY
        //----------------------------------------------------------------------

        printf ("\n") ;

        printf (         "Avg: NumberOfWalks all-pairs    k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_all [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: NumberOfWalks all-pairs    k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_all [nthreads], matrix_name) ;

        printf (         "Avg: MultiSourceBFS             sources: %" PRIu64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            (uint64_t) ntrials, nthreads, t_msbfs [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: MultiSourceBFS             sources: %" PRIu64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            (uint64_t) ntrials, nthreads, t_msbfs [nthreads], matrix_name) ;

        printf (         "Avg: NumberOfWalks single-src   k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_src [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: NumberOfWalks single-src   k: %2" PRId64
            "  threads: %3d  time: %10.3f sec  graph: %s\n",
            k, nthreads, t_walks_src [nthreads], matrix_name) ;

        printf (         "Avg: BFS single-source          threads: %3d"
            "  time: %10.3f sec  graph: %s\n",
            nthreads, t_bfs [nthreads], matrix_name) ;
        fprintf (stderr, "Avg: BFS single-source          threads: %3d"
            "  time: %10.3f sec  graph: %s\n",
            nthreads, t_bfs [nthreads], matrix_name) ;

        fflush (stdout) ; fflush (stderr) ;
    }

    // restore default thread count
    LAGRAPH_TRY (LAGraph_SetNumThreads (nthreads_outer, nthreads_inner, msg)) ;

    //--------------------------------------------------------------------------
    // CLEANUP
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
