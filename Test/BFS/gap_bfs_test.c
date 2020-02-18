//------------------------------------------------------------------------------
// bfs_test: read in (or create) a matrix and test BFS
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// usage:
// bfs_test s < in > out
// s is the staring node, in is the Matrix Market file, out is the level set.

// TODO clean this up

#include "../../Source/Utility/LAGraph_internal.h"
#include "bfs_test.h"

#define NTHREAD_LIST 6
#define THREAD_LIST 0

// #define NTHREAD_LIST 8
// #define THREAD_LIST 8, 7, 6, 5, 4, 3, 2, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&AT) ;                \
    GrB_free (&A) ;                 \
    GrB_free (&Abool) ;             \
    GrB_free (&v) ;                 \
    GrB_free (&pi) ;                \
    GrB_free (&SourceNodes) ;       \
}

int main (int argc, char **argv)
{
    GrB_Info info ;

    GrB_Matrix AT = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector v = NULL ;
    GrB_Vector pi = NULL ;
    GrB_Matrix SourceNodes = NULL ;
    LAGRAPH_OK (LAGraph_init ( )) ;

    uint64_t seed = 1 ;
    FILE *f ;
    int nthreads ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max = LAGraph_get_nthreads ( ) ;
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

    double t [nthreads_max+1] ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    double tic [2] ;
    LAGraph_tic (tic) ;

    double chunk ; // = 64 * 1024 ;
    GxB_get (GxB_CHUNK, &chunk) ;
    printf ("chunk: %g\n", chunk) ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    if (argc > 1)
    {
        // Usage:
        //      ./bfs_test matrixfile.mtx sources.mtx
        //      ./bfs_test matrixfile.grb sources.mtx

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }
        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            LAGRAPH_OK (LAGraph_binread (&A, filename)) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_mmread(&A, f));
            fclose (f) ;
        }

        // read in source nodes in Matrix Market format from the input file
        if (argc > 2)
        {
            filename = argv [2] ;
            printf ("sources: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Source node file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_mmread (&SourceNodes, f)) ;
            fclose (f) ;
        }
    }
    else
    {

        // Usage:  ./bfs_test < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_mmread(&A, stdin));
    }

    // convert to boolean, pattern-only
    LAGRAPH_OK (LAGraph_pattern (&Abool, A, GrB_BOOL)) ;
    // LAGraph_mmwrite (Abool, stderr) ;
    GrB_free (&A) ;
    A = Abool ;
    Abool = NULL ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols, nvals ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, A)) ;
    GrB_Index n = nrows ;
    if (nrows != ncols) { printf ("A must be square\n") ; abort ( ) ; }
    double t_read = LAGraph_toc (tic) ;
    printf ("read time: %g\n", t_read) ;

    // GxB_print (A, 2) ;

    //--------------------------------------------------------------------------
    // AT = A'
    //--------------------------------------------------------------------------

    // HACK: AT not needed for push-only BFS
    AT = NULL ;

#if 1
    LAGraph_tic (tic);
    bool A_is_symmetric ;
    LAGRAPH_OK (GrB_Matrix_new (&AT, GrB_BOOL, n, n)) ;
    LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
    LAGRAPH_OK (LAGraph_isequal (&A_is_symmetric, A, AT, NULL)) ;
    if (A_is_symmetric)
    {
        printf ("A is symmetric\n") ;
        GrB_free (&AT) ;
        AT = A ;
    }
    else
    {
        printf ("A is unsymmetric\n") ;
    }
    double t_transpose = LAGraph_toc (tic) ;
    printf ("transpose time: %g\n", t_transpose) ;
#endif

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    #define NSOURCES 64

    if (SourceNodes == NULL)
    {
        LAGRAPH_OK (GrB_Matrix_new (&SourceNodes, GrB_INT64, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            int64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
            // SourceNodes [k] = i 
            LAGRAPH_OK (GrB_Matrix_setElement (SourceNodes, i, k, 0)) ;
        }
    }

    int64_t ntrials ;
    GrB_Matrix_nrows (&ntrials, SourceNodes) ;

    printf ( "\n==========input graph: nodes: %lu edges: %lu ntrials: %lu\n",
        n, nvals, ntrials) ;

    //--------------------------------------------------------------------------
    // run the BFS on all source nodes
    //--------------------------------------------------------------------------

#if 0

    printf ( "simple all-push (no tree):\n") ;
    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_set_nthreads (nthreads) ;

        LAGraph_tic (tic) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (LAGraph_bfs_simple (&v, A, s)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:simple    (no tree): %12.3f (sec), rate: %6.2f\n",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (n > 1000)
        {
            LAGr_log (matrix_name, "simple:allpush", nthreads, t [nthreads]) ;
        }
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;

    char filename [1024] ;
    #if 0
    // dump the results so it can be checked
    sprintf (filename, "v_%d_simple.mtx", (int) n) ;
    f = fopen (filename, "w") ;
    LAGraph_mmwrite ((GrB_Matrix) v, f) ;
    fclose (f) ;
    #endif

    GrB_free (&v) ;

#endif

    //--------------------------------------------------------------------------
    // now the BFS on node s using push-pull (BEST) instead, no tree
    //--------------------------------------------------------------------------

#if 0
    printf ( "\n") ;
    for (int nthreads = 1 ; nthreads <= nthreads_max ; nthreads *= 2)
    {
        LAGraph_set_nthreads (nthreads) ;
        LAGraph_tic (tic) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v5) ;
            LAGRAPH_OK (LAGraph_bfs_pushpull (&v5, NULL, A, AT, s, 0, false)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:push/pull: (no tree) %12.3f (sec), "
            " rate: %6.2f\n",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;
    // make v5 sparse
    LAGRAPH_OK (GrB_assign (v5, v5, NULL, v5, GrB_ALL, n, LAGraph_desc_ooor)) ;
#endif

    //--------------------------------------------------------------------------
    // BFS: pushpull, with tree
    //--------------------------------------------------------------------------

    printf ( "pushpull (with tree):\n") ;
    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_set_nthreads (nthreads) ;
        t [nthreads] = 0 ;
        printf ("\n------------------------------------------- threads: %2d\n",
            nthreads) ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGraph_tic (tic) ;
            LAGRAPH_OK (LAGraph_bfs_pushpull (&v, &pi, A, AT, s, 0, false)) ;
            double ttrial = LAGraph_toc (tic) ;
            t [nthreads] += ttrial ;
            printf ("trial: %2d threads: %2d source: %9ld time: %10.4f sec\n",
                trial, nthreads, s, ttrial) ;
        }
        t [nthreads] = t [nthreads] / ntrials ;
        printf ( ":%2d:pushpull  (w/ tree): %12.3f (sec), rate: %6.2f\n",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (n > 1000)
        {
            LAGr_log (matrix_name, "w/tree:pushpull", nthreads, t [nthreads]) ;
        }
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;

    #if 0
    LAGraph_tic (tic) ;
    printf ("saving results ...\n")  ;

    // dump the last result so it can be checked
    sprintf (filename, "s_%d.mtx", (int) n) ;
    f = fopen (filename, "w") ;
    LAGraph_mmwrite ((GrB_Matrix) SourceNodes, f) ;
    fclose (f) ;

    sprintf (filename, "v_%d.mtx", (int) n) ;
    f = fopen (filename, "w") ;
    LAGraph_mmwrite ((GrB_Matrix) v, f) ;
    fclose (f) ;

    sprintf (filename, "p_%d.mtx", (int) n) ;
    f = fopen (filename, "w") ;
    LAGraph_mmwrite ((GrB_Matrix) pi, f) ;
    fclose (f) ;

    // LAGRAPH_OK (GxB_print (v, 2)) ;
    // LAGRAPH_OK (GxB_print (pi, 2)) ;
    GrB_free (&v) ;
    GrB_free (&pi) ;

    double t_save = LAGraph_toc (tic) ;
    printf ("save time: %g sec\n\n", t_save) ;
    #endif

    //--------------------------------------------------------------------------
    // BFS: both push and pull, with tree (log all timings)
    //--------------------------------------------------------------------------

#if 0
    printf ( "both (with tree, log all timings):\n") ;
    for (int tt = 1 ; tt <= nt ; tt++)
    {

        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_set_nthreads (nthreads) ;

        char logfilename [1024] ;
        sprintf (logfilename, "logtime_%ld_%d.txt", n, nthreads) ;
        printf ("logfile: [%s]\n", logfilename) ;
        FILE *logfile = fopen (logfilename, "w") ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (LAGraph_bfs_both (&v, &pi, A, AT, s, 0, false,
                logfile)) ;
        }

        fclose (logfile) ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;
#endif

#if 0

    sprintf (filename, "allpush_%lu.m", n) ;
    f = fopen (filename, "w") ;
    fprintf (f, "function [results, name] = allpush_%lu\n", n);
    fprintf (f, "%% bfs_log: all-push\n") ;
    fprintf (f, "name = '%s' ;\n", matrix_name) ;
    fprintf (f, "n = %lu ;\n", n) ;
    fprintf (f, "nvals = %lu ;\n", nvals) ;
    fprintf (f, "d = %g ;\n", ((double) nvals) / (double) n) ;
    fprintf (f, "%%%% columns in results:\n") ;
    fprintf (f, "%%%% s, n, d, nthreads, do_push, level, nq, nvisited, time\n");
    fprintf (f, "results = [\n") ;
    printf ( "allpush (log, with tree):\n") ;
    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_set_nthreads (nthreads) ;
        LAGraph_tic (tic) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (bfs_log (&v, &pi, A, NULL, s, 0, false, f)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:allpush   (w/ tree): %12.3f (sec), rate: %6.2f\n",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;
    // TODO: check results
    // LAGRAPH_OK (GxB_print (v, 2)) ;
    // LAGRAPH_OK (GxB_print (pi, 2)) ;
    GrB_free (&v) ;
    GrB_free (&pi) ;
    fprintf (f, "] ;\n") ;
    fclose (f) ;

#endif

    //--------------------------------------------------------------------------
    // BFS: all-pull, with tree (log all timings)
    //--------------------------------------------------------------------------

#if 0

    nthreads = nthreads_max ;
    sprintf (filename, "allpull_%lu.m", n) ;
    f = fopen (filename, "w") ;
    fprintf (f, "function [results, name] = allpush_%lu\n", n);
    fprintf (f, "%% bfs_log: all-pull\n") ;
    fprintf (f, "name = '%s' ;\n", matrix_name) ;
    fprintf (f, "n = %lu ;\n", n) ;
    fprintf (f, "nvals = %lu ;\n", nvals) ;
    fprintf (f, "d = %g ;\n", ((double) nvals) / (double) n) ;
    fprintf (f, "%%%% columns in results:\n") ;
    fprintf (f, "%%%% s, n, d, nthreads, do_push, level, nq, nvisited, time\n");
    fprintf (f, "results = [\n") ;
    printf ( "allpull (log, with tree):\n") ;
    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_set_nthreads (nthreads) ;
        LAGraph_tic (tic) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (bfs_log (&v, &pi, NULL, AT, s, 0, false, f)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:allpull   (w/ tree): %12.3f (sec), rate: %6.2f\n",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;
    // TODO: check results
    // LAGRAPH_OK (GxB_print (v, 2)) ;
    // LAGRAPH_OK (GxB_print (pi, 2)) ;
    GrB_free (&v) ;
    GrB_free (&pi) ;
    fprintf (f, "] ;\n") ;
    fclose (f) ;

#endif

    //--------------------------------------------------------------------------
    // BFS: push-pull, with tree
    //--------------------------------------------------------------------------

#if 0

    printf ( "pushpull (with tree):\n") ;
    for (int tt = 1 ; tt <= nt ; tt++)
    {
        int nthreads = Nthreads [tt] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_set_nthreads (nthreads) ;
        LAGraph_tic (tic) ;
        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t s ; 
            // s = SourceNodes [i]
            LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
            s-- ; // convert from 1-based to 0-based
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (LAGraph_bfs_pushpull (&v, &pi, A, AT, s, 0, false)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:push/pull (w/ tree): %12.3f (sec), rate: %6.2f\n",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;
    // TODO: check results
    // LAGRAPH_OK (GxB_print (v, 2)) ;
    // LAGRAPH_OK (GxB_print (pi, 2)) ;
    GrB_free (&v) ;
    GrB_free (&pi) ;

#endif

#if 0

    //--------------------------------------------------------------------------
    // now the BFS on node s using push-pull with v sparse
    //--------------------------------------------------------------------------

    printf ( "v starts sparse: (nthread %d)\n", nthreads_max) ;
    LAGraph_tic (tic) ;
    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        int64_t s ; 
        // s = SourceNodes [i]
        LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
        s-- ; // convert from 1-based to 0-based
        GrB_free (&v6) ;
        LAGRAPH_OK (LAGraph_bfs_pushpull (&v6, NULL, A, AT, s, 0, true)) ;
    }
    double t6 = LAGraph_toc (tic) / ntrials ;
    printf ( "pushpull sparse %12.3f (sec), rate: %6.2f\n",
        t6, 1e-6*((double) nvals) / t6) ;

    int32_t maxlevel = 0 ;
#endif

    //--------------------------------------------------------------------------
    // check results
    //--------------------------------------------------------------------------

#if 0
    // to make result sparse
    // LAGRAPH_OK (GrB_assign (v, v, NULL, v, GrB_ALL, n, LAGraph_desc_ooor)) ;

    bool isequal = false, ok = true ;
    int64_t s ; 
    LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, ntrials-1, 0)) ;
    s-- ; // convert from 1-based to 0-based

    // find the max level
    LAGRAPH_OK (GrB_reduce (&maxlevel, NULL, LAGraph_MAX_INT32_MONOID, v,
        NULL));
    printf ( "number of levels: %d (for s = %lu, single-source)\n",
        maxlevel, s) ;

    // find the number of nodes visited
    GrB_Index nv = 0 ;
    LAGRAPH_OK (GrB_Vector_nvals (&nv, v)) ;
    printf ( "# nodes visited (for single-source): %lu out of %lu"
        " (%g %% of the graph)\n", nv, n,
        100. * (double) nv / (double) n) ;

    LAGRAPH_OK (LAGraph_Vector_isequal (&isequal, v, v5_tree, NULL)) ;
    if (!isequal)
    {
        printf ( "ERROR! simple and TREE   differ\n") ;
        ok = false ;
    }

    LAGRAPH_OK (LAGraph_Vector_isequal (&isequal, v, v5, NULL)) ;
    if (!isequal)
    {
        printf ( "ERROR! simple and best   differ\n") ;
        ok = false ;
    }

    LAGRAPH_OK (LAGraph_Vector_isequal (&isequal, v, v6, NULL)) ;
    if (!isequal)
    {
        // GxB_print (v, 2) ;
        // GxB_print (v6, 2) ;
        printf ( "ERROR! simple and push-pull (sparse) differ\n") ;
        ok = false ;
    }
#endif

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_finalize ( )) ;

#if 0
    printf ( "bfs_test: ") ;
    if (ok)
    {
        printf ( "all tests passed\n") ;
    }
    else
    {
        printf ( "TEST FAILURE\n") ;
    }
    printf ("----------------------------------------------------------\n\n") ;
#endif

    return (GrB_SUCCESS) ;
}

