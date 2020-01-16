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
#define THREAD_LIST 64, 32, 24, 12, 8, 4

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
    uint64_t seed = 1 ;
    FILE *f ;

    GrB_Matrix AT = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector v = NULL ;
    GrB_Vector pi = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    LAGRAPH_OK (LAGraph_init ( )) ;
    int nthreads ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max = LAGraph_get_nthreads();
    Nthreads [nt] = LAGRAPH_MIN (Nthreads [nt], nthreads_max) ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" thread test %d: %d\n", t, nthreads) ;
    }

    double t [nthreads_max+1] ;
    char filename [1000] ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    double tic [2] ;
    LAGraph_tic (tic) ;

    double chunk ; // = 64 * 1024 ;
    GxB_get (GxB_CHUNK, chunk) ;
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

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    if (SourceNodes == NULL)
    {
        LAGRAPH_OK (GrB_Matrix_new (&SourceNodes, GrB_INT64, 64, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < 64 ; k++)
        {
            int64_t i = rand ( ) % n ;
            // SourceNodes [k] = i 
            LAGRAPH_OK (GrB_Matrix_setElement (SourceNodes, i, k, 0)) ;
        }
    }

    // GxB_print (SourceNodes, 3) ;

    int64_t ntrials ;
    GrB_Matrix_nrows (&ntrials, SourceNodes) ;

    printf ( "\n==========input graph: nodes: %lu edges: %lu ntrials: %lu\n",
        n, nvals, ntrials) ;

    //--------------------------------------------------------------------------
    // run the BFS on all source nodes
    //--------------------------------------------------------------------------

#if 0
    int nth2 = LAGraph_get_nthreads ( ) ;
    printf ("nthreads for simple BFS: %d\n", nth2) ;

    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        int64_t s ; 
        // s = SourceNodes [i]
        LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0)) ;
        GrB_free (&v) ;
        LAGRAPH_OK (LAGraph_bfs_simple (&v, A, s)) ;
    }

    // stop the timer
    double t1 = LAGraph_toc (tic) / ntrials ;
    printf ( "nthreads_max: %d\n", nthreads_max) ;
    printf ( "simple    time: %12.3f (sec), rate: %6.2f\n",
        t1, 1e-6*((double) nvals) / t1) ;
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
            GrB_free (&v5) ;
            LAGRAPH_OK (LAGraph_bfs_pushpull (&v5, NULL, A, AT, s, 0, false)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:push/pull: (no tree) %12.3f (sec), "
            " rate: %6.2f",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (nthreads > 1) printf ( " speedup %4.1f", t [1] / t [nthreads]) ;
        printf ("\n") ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;
    // make v5 sparse
    LAGRAPH_OK (GrB_assign (v5, v5, NULL, v5, GrB_ALL, n, LAGraph_desc_ooor)) ;
#endif

    //--------------------------------------------------------------------------
    // BFS: all-push, with tree
    //--------------------------------------------------------------------------

    printf ( "allpush (with tree):\n") ;
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
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (LAGraph_bfs_pushpull (&v, &pi, A, NULL, s, 0, false)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:allpush   (w/ tree): %12.3f (sec), rate: %6.2f",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (nthreads > 1) printf ( " speedup %4.1f", t [1] / t [nthreads]) ;
        printf ("\n") ;
    }
    // restore default
    LAGraph_set_nthreads (nthreads_max) ;
    printf ( "\n") ;

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

    //--------------------------------------------------------------------------
    // BFS: all-push, with tree (log all timings)
    //--------------------------------------------------------------------------

#if 1

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
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (bfs_log (&v, &pi, A, NULL, s, 0, false, f)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:allpush   (w/ tree): %12.3f (sec), rate: %6.2f",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (nthreads > 1) printf ( " speedup %4.1f", t [1] / t [nthreads]) ;
        printf ("\n") ;
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

#if 1

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
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (bfs_log (&v, &pi, NULL, AT, s, 0, false, f)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:allpull   (w/ tree): %12.3f (sec), rate: %6.2f",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (nthreads > 1) printf ( " speedup %4.1f", t [1] / t [nthreads]) ;
        printf ("\n") ;
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
            GrB_free (&v) ;
            GrB_free (&pi) ;
            LAGRAPH_OK (LAGraph_bfs_pushpull (&v, &pi, A, AT, s, 0, false)) ;
        }
        t [nthreads] = LAGraph_toc (tic) / ntrials ;
        printf ( ":%2d:push/pull (w/ tree): %12.3f (sec), rate: %6.2f",
            nthreads, t [nthreads], 1e-6*((double) nvals) / t [nthreads]) ;
        if (nthreads > 1) printf ( " speedup %4.1f", t [1] / t [nthreads]) ;
        printf ("\n") ;
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

