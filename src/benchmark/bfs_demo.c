//------------------------------------------------------------------------------
// LAGraph/src/benchmark/bfs_demo.c: benchmark for LAGr_BreadthFirstSearch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LAGraph_demo.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0

// #define NTHREAD_LIST 8
// #define THREAD_LIST 8, 7, 6, 5, 4, 3, 2, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, msg) ;      \
    GrB_free (&A) ;                 \
    GrB_free (&Abool) ;             \
    GrB_free (&parent) ;            \
    GrB_free (&level) ;             \
    GrB_free (&SourceNodes) ;       \
}

int main (int argc, char **argv)
{

    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector level = NULL ;
    GrB_Vector parent = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    uint64_t seed = 1 ;
    FILE *f ;
    int nthreads ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
    printf ("nthreads_max: %d\n", nthreads_max) ;
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

    double tpl [nthreads_max+1] ;
    double tp [nthreads_max+1] ;
    double tl [nthreads_max+1] ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, &SourceNodes,
        false, false, true, NULL, false, argc, argv)) ;

    // compute G->rowdegree
    LAGRAPH_TRY (LAGraph_Property_RowDegree (G, msg)) ;

    // compute G->coldegree, just to test it (not needed for any tests)
    LAGRAPH_TRY (LAGraph_Property_ColDegree (G, msg)) ;

    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    GrB_Index ntrials ;
    GRB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;

    // HACK
    // ntrials = 4 ;

    //--------------------------------------------------------------------------
    // warmup
    //--------------------------------------------------------------------------

    int64_t src ;
    double twarmup, tw [2] ;
    GRB_TRY (GrB_Matrix_extractElement (&src, SourceNodes, 0, 0)) ;
    LAGRAPH_TRY (LAGraph_Tic (tw, msg)) ;
    LAGRAPH_TRY (LAGr_BreadthFirstSearch (NULL, &parent, G, src, msg)) ;
    GrB_free (&parent) ;
    LAGRAPH_TRY (LAGraph_Toc (&twarmup, tw, msg)) ;
    printf ("warmup: parent only, pushpull: %g sec\n", twarmup) ;

    //--------------------------------------------------------------------------
    // run the BFS on all source nodes
    //--------------------------------------------------------------------------

    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (nthreads, msg)) ;

        tp [nthreads] = 0 ;
        tl [nthreads] = 0 ;
        tpl [nthreads] = 0 ;

        printf ("\n------------------------------- threads: %2d\n", nthreads) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t src ;
            // src = SourceNodes [trial]
            GRB_TRY (GrB_Matrix_extractElement (&src, SourceNodes, trial, 0)) ;
            src-- ; // convert from 1-based to 0-based
            double tcheck, ttrial, tic [2] ;

            {

                //--------------------------------------------------------------
                // BFS to compute just parent
                //--------------------------------------------------------------

                GrB_free (&parent) ;
                LAGRAPH_TRY (LAGraph_Tic (tic, msg)) ;
                LAGRAPH_TRY (LAGr_BreadthFirstSearch (NULL, &parent,
                    G, src, msg)) ;
                LAGRAPH_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
                tp [nthreads] += ttrial ;
                printf ("parent only  pushpull trial: %2d threads: %2d "
                    "src: %g %10.4f sec\n",
                    trial, nthreads, (double) src, ttrial) ;
                fflush (stdout) ;

                int32_t maxlevel ;
                GrB_Index nvisited ;

#if LG_CHECK_RESULT
                // check the result (this is very slow so only do it for one trial)
                if (trial == 0)
                {
                    LAGRAPH_TRY (LAGraph_Tic (tic, msg)) ;
                    LAGRAPH_TRY (LG_check_bfs (NULL, parent, G, src, msg)) ;
                    LAGRAPH_TRY (LAGraph_Toc (&tcheck, tic, msg)) ;
                    printf ("    n: %g check: %g sec\n", (double) n, tcheck) ;
                }
#endif

                GrB_free (&parent) ;

                //--------------------------------------------------------------
                // BFS to compute just level
                //--------------------------------------------------------------

#if 0
                GrB_free (&level) ;

                LAGRAPH_TRY (LAGraph_Tic (tic, msg)) ;
                LAGRAPH_TRY (LAGr_BreadthFirstSearch (&level, NULL,
                    G, src, msg)) ;
                LAGRAPH_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
                tl [nthreads] += ttrial ;

                GRB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32,
                    level, NULL)) ;
                printf ("level only   pushpull trial: %2d threads: %2d "
                    "src: %g %10.4f sec maxlevel: %d\n",
                    trial, nthreads, (double) src, ttrial, maxlevel) ;
                fflush (stdout) ;

#if LG_CHECK_RESULT
                // check the result (this is very slow so only do it for one trial)
                if (trial == 0)
                {
                    LAGRAPH_TRY (LAGraph_Tic (tic, msg)) ;
                    LAGRAPH_TRY (LG_check_bfs (level, NULL, G, src, msg)) ;
                    GRB_TRY (GrB_Vector_nvals (&nvisited, level)) ;
                    LAGRAPH_TRY (LAGraph_Toc (&tcheck, tic, msg)) ;
                    printf ("    n: %g max level: %d nvisited: %g "
                        "check: %g sec\n", (double) n, maxlevel,
                        (double) nvisited, tcheck) ;
                }
#endif

                GrB_free (&level) ;

                //--------------------------------------------------------------
                // BFS to compute both parent and level
                //--------------------------------------------------------------

                GrB_free (&parent) ;
                GrB_free (&level) ;
                LAGRAPH_TRY (LAGraph_Tic (tic, msg)) ;
                LAGRAPH_TRY (LAGr_BreadthFirstSearch (&level, &parent,
                    G, src, msg)) ;
                LAGRAPH_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
                tpl [nthreads] += ttrial ;

                GRB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32,
                    level, NULL)) ;
                printf ("parent+level pushpull trial: %2d threads: %2d "
                    "src: %g %10.4f sec\n",
                    trial, nthreads, (double) src, ttrial) ;
                fflush (stdout) ;

#if LG_CHECK_RESULT
                // check the result (this is very slow so only do it for one trial)
                if (trial == 0)
                {
                    LAGRAPH_TRY (LAGraph_Tic (tic, msg)) ;
                    LAGRAPH_TRY (LG_check_bfs (level, parent, G, src, msg)) ;
                    GRB_TRY (GrB_Vector_nvals (&nvisited, level)) ;
                    LAGRAPH_TRY (LAGraph_Toc (&tcheck, tic, msg)) ;
                    printf ("    n: %g max level: %d nvisited: %g "
                        "check: %g sec\n",
                        (double) n, maxlevel, (double) nvisited, tcheck) ;
                }
#endif
#endif

                GrB_free (&parent) ;
                GrB_free (&level) ;
            }
        }

        {
            tp  [nthreads] = tp  [nthreads] / ntrials ;
            tl  [nthreads] = tl  [nthreads] / ntrials ;
            tpl [nthreads] = tpl [nthreads] / ntrials ;

            fprintf (stderr, "Avg: BFS pushpull parent only  threads %3d: "
                "%10.3f sec: %s\n",
                 nthreads, tp [nthreads], matrix_name) ;
#if 0
            fprintf (stderr, "Avg: BFS pushpull level only   threads %3d: "
                "%10.3f sec: %s\n",
                 nthreads, tl [nthreads], matrix_name) ;

            fprintf (stderr, "Avg: BFS pushpull level+parent threads %3d: "
                "%10.3f sec: %s\n",
                 nthreads, tpl [nthreads], matrix_name) ;
#endif

            printf ("Avg: BFS pushpull parent only  threads %3d: "
                "%10.3f sec: %s\n",
                 nthreads, tp [nthreads], matrix_name) ;

#if 0
            printf ("Avg: BFS pushpull level only   threads %3d: "
                "%10.3f sec: %s\n",
                 nthreads, tl [nthreads], matrix_name) ;

            printf ("Avg: BFS pushpull level+parent threads %3d: "
                "%10.3f sec: %s\n",
                 nthreads, tpl [nthreads], matrix_name) ;
#endif
        }
    }
    // restore default
    LAGRAPH_TRY (LAGraph_SetNumThreads (nthreads_max, msg)) ;
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
