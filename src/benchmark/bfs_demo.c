//------------------------------------------------------------------------------
// LAGraph/Test2/BreadthFirstSearch/test_bfs.c: test LAGraph_BreadthFirstSearch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_demo.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0

// #define NTHREAD_LIST 8
// #define THREAD_LIST 8, 7, 6, 5, 4, 3, 2, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGraph_FREE_ALL            \
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
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
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

    double tpl [nthreads_max+1][2] ;
    double tp [nthreads_max+1][2] ;
    double tl [nthreads_max+1][2] ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    if (readproblem (&G, &SourceNodes,
        false, false, true, NULL, false, argc, argv) != 0) ERROR ;

    // compute G->rowdegree
    LAGraph_TRY (LAGraph_Property_RowDegree (G, msg)) ;

    // compute G->coldegree, just to test it (not needed for any tests)
    LAGraph_TRY (LAGraph_Property_ColDegree (G, msg)) ;

    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    int64_t ntrials ;
    GrB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;

    // HACK
    // ntrials = 4 ;

    //--------------------------------------------------------------------------
    // run the BFS on all source nodes
    //--------------------------------------------------------------------------

    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_TRY (LAGraph_SetNumThreads (nthreads, msg)) ;

        tp [nthreads][0] = 0 ;
        tl [nthreads][0] = 0 ;
        tpl [nthreads][0] = 0 ;

        tp [nthreads][1] = 0 ;
        tl [nthreads][1] = 0 ;
        tpl [nthreads][1] = 0 ;

        printf ("\n------------------------------- threads: %2d\n", nthreads) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t src ;
            // src = SourceNodes [trial]
            GrB_TRY (GrB_Matrix_extractElement (&src, SourceNodes, trial, 0)) ;
            src-- ; // convert from 1-based to 0-based
            double tcheck, ttrial, tic [2] ;

            for (int pp = 0 ; pp <= 1 ; pp++)
            {

                bool pushpull = (pp == 1) ;

                //--------------------------------------------------------------
                // BFS to compute just parent
                //--------------------------------------------------------------

                GrB_free (&parent) ;
                LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
                LAGraph_TRY (LAGraph_BreadthFirstSearch (NULL, &parent,
                    G, src, pushpull, msg)) ;
                LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
                tp [nthreads][pp] += ttrial ;
                printf ("parent only  %s trial: %2d threads: %2d "
                    "src: %9ld %10.4f sec\n",
                    (pp == 0) ? "pushonly" : "pushpull",
                    trial, nthreads, src, ttrial) ;
                fflush (stdout) ;

                // check the result (this is very slow so only do it for one trial)
                int32_t maxlevel ;
                int64_t nvisited ;
                if (trial == 0)
                {
                    LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
                    LAGraph_TRY (LG_check_bfs (NULL, parent, G, src, msg)) ;
                    LAGraph_TRY (LAGraph_Toc (&tcheck, tic, msg)) ;
                    printf ("    n: %ld check: %g sec\n", n, tcheck) ;
                }

                GrB_free (&parent) ;

                //--------------------------------------------------------------
                // BFS to compute just level
                //--------------------------------------------------------------

#if 1
                GrB_free (&level) ;

                LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
                LAGraph_TRY (LAGraph_BreadthFirstSearch (&level, NULL,
                    G, src, pushpull, msg)) ;
                LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
                tl [nthreads][pp] += ttrial ;

                GrB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32,
                    level, NULL)) ;
                printf ("level only   %s trial: %2d threads: %2d "
                    "src: %9ld %10.4f sec\n",
                    (pp == 0) ? "pushonly" : "pushpull",
                    trial, nthreads, src, ttrial) ;
                fflush (stdout) ;

                // check the result (this is very slow so only do it for one trial)
                if (trial == 0)
                {
                    LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
                    LAGraph_TRY (LG_check_bfs (level, NULL, G, src, msg)) ;
                    GrB_TRY (GrB_Vector_nvals (&nvisited, level)) ;
                    LAGraph_TRY (LAGraph_Toc (&tcheck, tic, msg)) ;
                    printf ("    n: %ld max level: %d nvisited: %ld check: %g sec\n",
                        n, maxlevel, nvisited, tcheck) ;
                }

                GrB_free (&level) ;

                //--------------------------------------------------------------
                // BFS to compute both parent and level
                //--------------------------------------------------------------

                GrB_free (&parent) ;
                GrB_free (&level) ;
                LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
                LAGraph_TRY (LAGraph_BreadthFirstSearch (&level, &parent,
                    G, src, pushpull, msg)) ;
                LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
                tpl [nthreads][pp] += ttrial ;

                GrB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32,
                    level, NULL)) ;
                printf ("parent+level %s trial: %2d threads: %2d "
                    "src: %9ld %10.4f sec\n",
                    (pp == 0) ? "pushonly" : "pushpull",
                    trial, nthreads, src, ttrial) ;
                fflush (stdout) ;

                // check the result (this is very slow so only do it for one trial)
                if (trial == 0)
                {
                    LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
                    LAGraph_TRY (LG_check_bfs (level, parent, G, src, msg)) ;
                    GrB_TRY (GrB_Vector_nvals (&nvisited, level)) ;
                    LAGraph_TRY (LAGraph_Toc (&tcheck, tic, msg)) ;
                    printf ("    n: %ld max level: %d nvisited: %ld check: %g sec\n",
                        n, maxlevel, nvisited, tcheck) ;
                }

#endif

                GrB_free (&parent) ;
                GrB_free (&level) ;
            }
        }

        for (int pp = 0 ; pp <= 1 ; pp++)
        {
            tp  [nthreads][pp] = tp  [nthreads][pp] / ntrials ;
            tl  [nthreads][pp] = tl  [nthreads][pp] / ntrials ;
            tpl [nthreads][pp] = tpl [nthreads][pp] / ntrials ;

            fprintf (stderr, "Avg: BFS %s parent only  threads %3d: "
                "%10.3f sec: %s\n",
                 (pp == 0) ? "pushonly" : "pushpull",
                 nthreads, tp [nthreads][pp], matrix_name) ;
#if 1
            fprintf (stderr, "Avg: BFS %s level only   threads %3d: "
                "%10.3f sec: %s\n",
                 (pp == 0) ? "pushonly" : "pushpull",
                 nthreads, tl [nthreads][pp], matrix_name) ;

            fprintf (stderr, "Avg: BFS %s level+parent threads %3d: "
                "%10.3f sec: %s\n",
                 (pp == 0) ? "pushonly" : "pushpull",
                 nthreads, tpl [nthreads][pp], matrix_name) ;
#endif

            printf ("Avg: BFS %s parent only  threads %3d: "
                "%10.3f sec: %s\n",
                 (pp == 0) ? "pushonly" : "pushpull",
                 nthreads, tp [nthreads][pp], matrix_name) ;

#if 1
            printf ("Avg: BFS %s level only   threads %3d: "
                "%10.3f sec: %s\n",
                 (pp == 0) ? "pushonly" : "pushpull",
                 nthreads, tl [nthreads][pp], matrix_name) ;

            printf ("Avg: BFS %s level+parent threads %3d: "
                "%10.3f sec: %s\n",
                 (pp == 0) ? "pushonly" : "pushpull",
                 nthreads, tpl [nthreads][pp], matrix_name) ;
#endif
        }
    }
    // restore default
    LAGraph_TRY (LAGraph_SetNumThreads (nthreads_max, msg)) ;
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGraph_FREE_ALL ;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}
