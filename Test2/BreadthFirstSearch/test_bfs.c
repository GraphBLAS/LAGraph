//------------------------------------------------------------------------------
// LAGraph/Test2/BreadthFirstSearch/test_bfs.c: test LAGraph_BreadthFirstSearch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_Test.h"

// #define NTHREAD_LIST 2
#define NTHREAD_LIST 1
#define THREAD_LIST 0

// #define NTHREAD_LIST 8
// #define THREAD_LIST 8, 7, 6, 5, 4, 3, 2, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGRAPH_FREE_ALL            \
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

    printf ("%s v%d.%d.%d [%s]\n",
        GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR,
        GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB,
        GxB_IMPLEMENTATION_DATE) ;

    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector level = NULL ;
    GrB_Vector parent = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    LAGraph_TRY (LAGraph_Init (msg)) ;
    GrB_TRY (GxB_set (GxB_BURBLE, false)) ;

    uint64_t seed = 1 ;
    FILE *f ;
    int nthreads ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ; 
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
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
    LAGraph_TRY (LAGraph_Test_ReadProblem (&G, &SourceNodes,
        false, false, true, NULL, false, argc, argv, msg)) ;

    // compute G->rowdegree
    LAGraph_TRY (LAGraph_Property_RowDegree (G, msg)) ;

    // compute G->coldegree, just to test it (not needed for any tests)
    LAGraph_TRY (LAGraph_Property_ColDegree (G, msg)) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    int64_t ntrials ;
    GrB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;

    // HACK
    // ntrials = 1 ;

    //--------------------------------------------------------------------------
    // run the BFS on all source nodes
    //--------------------------------------------------------------------------

    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_TRY (LAGraph_SetNumThreads (nthreads, msg)) ;
        tp [nthreads] = 0 ;
        tl [nthreads] = 0 ;
        tpl [nthreads] = 0 ;
        printf ("\n------------------------------- threads: %2d\n", nthreads) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t src ; 
            // src = SourceNodes [trial]
            GrB_TRY (GrB_Matrix_extractElement (&src, SourceNodes, trial, 0)) ;
            src-- ; // convert from 1-based to 0-based
            double ttrial, tic [2] ;

            //------------------------------------------------------------------
            // BFS to compute just parent
            //------------------------------------------------------------------

            GrB_free (&parent) ;
            LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
            LAGraph_TRY (LAGraph_BreadthFirstSearch (NULL, &parent,
                G, src, msg)) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
            tp [nthreads] += ttrial ;
            printf ("parent only  trial: %2d threads: %2d src: %9ld "
                "%10.4f sec\n", trial, nthreads, src, ttrial) ;
            fflush (stdout) ;
            // GrB_TRY (GxB_print (parent, 2)) ;
            GrB_free (&parent) ;

            //------------------------------------------------------------------
            // BFS to compute just level
            //------------------------------------------------------------------

            GrB_free (&level) ;
#if 1
            LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
            LAGraph_TRY (LAGraph_BreadthFirstSearch (&level, NULL,
                G, src, msg)) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
            tl [nthreads] += ttrial ;

            int32_t maxlevel ;
            GrB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32, level,
                NULL)) ;

            printf ("level only   trial: %2d threads: %2d src: %9ld "
                "%10.4f sec maxlevel %d\n",
                trial, nthreads, src, ttrial, maxlevel) ;
            fflush (stdout) ;
            // GrB_TRY (GxB_print (level, 2)) ;

            GrB_free (&level) ;

            //------------------------------------------------------------------
            // BFS to compute both parent and level
            //------------------------------------------------------------------

            GrB_free (&parent) ;
            GrB_free (&level) ;
            LAGraph_TRY (LAGraph_Tic (tic, msg)) ;
            LAGraph_TRY (LAGraph_BreadthFirstSearch (&level, &parent,
                G, src, msg)) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, msg)) ;
            tpl [nthreads] += ttrial ;

            GrB_TRY (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT32, level,
                NULL)) ;
            printf ("parent+level trial: %2d threads: %2d src: %9ld "
                "%10.4f sec maxlevel %d\n",
                trial, nthreads, src, ttrial, maxlevel) ;
            fflush (stdout) ;
            // GrB_TRY (GxB_print (parent, 2)) ;
            // GrB_TRY (GxB_print (level, 2)) ;
#endif
            GrB_free (&parent) ;
            GrB_free (&level) ;
        }

        tp [nthreads] = tp [nthreads] / ntrials ;
        tl [nthreads] = tl [nthreads] / ntrials ;
        tpl [nthreads] = tpl [nthreads] / ntrials ;

        fprintf (stderr, "Avg: BFS parent only  %3d: %10.3f sec: %s\n",
             nthreads, tp [nthreads], matrix_name) ;

        fprintf (stderr, "Avg: BFS level only   %3d: %10.3f sec: %s\n",
             nthreads, tl [nthreads], matrix_name) ;

        fprintf (stderr, "Avg: BFS level+parent %3d: %10.3f sec: %s\n",
             nthreads, tpl [nthreads], matrix_name) ;

    }
    // restore default
    LAGraph_TRY (LAGraph_SetNumThreads (nthreads_max, msg)) ;
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}

