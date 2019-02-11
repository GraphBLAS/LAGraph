//------------------------------------------------------------------------------
// bfs_test: read in (or create) a matrix and test BFS
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// usage:
// bfs_test s < in > out
// s is the staring node, in is the Matrix Market file, out is the level set.

// TODO clean this up

#include "bfs_test.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&AT) ;                \
    GrB_free (&A) ;                 \
    GrB_free (&Abool) ;             \
    GrB_free (&v) ;                 \
    GrB_free (&v2) ;                \
    GrB_free (&v3) ;                \
    GrB_free (&v4) ;                \
    GrB_free (&v5) ;                \
    GrB_free (&diff) ;              \
}

int main (int argc, char **argv)
{

    GrB_Info info ;
    uint64_t seed = 1 ;

    GrB_Matrix AT = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector v = NULL ;
    GrB_Vector v2 = NULL ;
    GrB_Vector v3 = NULL ;
    GrB_Vector v4 = NULL ;
    GrB_Vector v5 = NULL ;
    GrB_Vector diff = NULL ;

    LAGRAPH_OK (LAGraph_init ( )) ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // read in the file in Matrix Market format
    LAGRAPH_OK (LAGraph_mmread (&A, stdin)) ;
    // GxB_fprint (A, GxB_COMPLETE, stderr) ;
    // LAGraph_mmwrite (A, stderr) ;

    // convert to boolean, pattern-only
    LAGRAPH_OK (LAGraph_pattern (&Abool, A)) ;
    // LAGraph_mmwrite (Abool, stderr) ;
    GrB_free (&A) ;
    A = Abool ;
    // GxB_fprint (A, GxB_COMPLETE, stderr) ;

    // finish any pending computations
    GrB_Index nvals ;
    GrB_Matrix_nvals (&nvals, A) ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;

    // GxB_fprint (A, GxB_COMPLETE, stderr) ;

    // LAGRAPH_OK (GrB_Matrix_setElement (A, 0, 0, n-1)) ;     // hack

    //--------------------------------------------------------------------------
    // get the source node
    //--------------------------------------------------------------------------

    GrB_Index s = 0 ;
    if (argc > 1)
    {
        sscanf (argv [1], "%" PRIu64, &s) ; 
    }

    fprintf (stderr, "\n=========="
        "input graph: nodes: %g edges: %g source node: %g\n",
        (double) n, (double)  nvals, (double) s) ;

    //--------------------------------------------------------------------------
    // run the BFS on node s
    //--------------------------------------------------------------------------

    int ntrials = 1000 ;

    // start the timer
    double tic [2] ;
    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&v) ;
        LAGRAPH_OK (LAGraph_bfs_simple (&v, A, s)) ;
    }

    // stop the timer
    double t1 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "simple    time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t1, 1e-6*((double) nvals) / t1) ;

    //--------------------------------------------------------------------------
    // run the BFS on node s with LAGraph_bfs2 (PUSH)
    //--------------------------------------------------------------------------

    // start the timer
    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&v2) ;
        LAGRAPH_OK (LAGraph_bfs2 (&v2, A, s, INT32_MAX)) ;
    }

    // stop the timer
    double t2 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "method2   time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t2, 1e-6*((double) nvals) / t2) ;
    fprintf (stderr, "speedup of method2:   %g\n", t1/t2) ;

    //--------------------------------------------------------------------------
    // AT = A'
    //--------------------------------------------------------------------------

    LAGraph_tic (tic) ;
    LAGRAPH_OK (GrB_Matrix_new (&AT, GrB_BOOL, ncols, nrows)) ;
    LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
    double transpose_time = LAGraph_toc (tic) ;
    fprintf (stderr, "transpose time: %g\n", transpose_time) ;

    //--------------------------------------------------------------------------
    // now the BFS on node s using push-pull instead
    //--------------------------------------------------------------------------

    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&v3) ;
        LAGRAPH_OK (LAGraph_bfs_pushpull_old (&v3, A, AT, s, INT32_MAX)) ;
    }

    double t3 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "push/pull old:  %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t3, 1e-6*((double) nvals) / t3) ;
    fprintf (stderr, "speedup of push/pull OLD: %g\n", t1/t3) ;

    //--------------------------------------------------------------------------
    // now the BFS on node s using push-pull (BEST) instead
    //--------------------------------------------------------------------------

    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&v5) ;
        LAGRAPH_OK (LAGraph_bfs_pushpull (&v5, A, AT, s, 0, false)) ;
    }

    double t5 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "push/pull best: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t5, 1e-6*((double) nvals) / t5) ;
    fprintf (stderr, "speedup of push/pull (best): %g\n", t1/t5) ;

    //--------------------------------------------------------------------------
    // now the BFS on node s using PULL (only) instead
    //--------------------------------------------------------------------------

    // slow!!!
    LAGraph_tic (tic) ;
    LAGRAPH_OK (LAGraph_bfs_pull (&v4, A, AT, s, INT32_MAX)) ;
    double t4 = LAGraph_toc (tic) ;
    fprintf (stderr, "pull      time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t4, 1e-6*((double) nvals) / t4) ;
    fprintf (stderr, "speedup of push/pull: %g (normally slow!)\n", t1/t4) ;

    //--------------------------------------------------------------------------
    // check results
    //--------------------------------------------------------------------------

    // find the max level
    int32_t maxlevel = 0 ;
    LAGRAPH_OK (GrB_reduce (&maxlevel, NULL, LAGraph_MAX_INT32_MONOID, v,
        NULL));
    fprintf (stderr, "number of levels: %d\n", maxlevel) ;

    // find the number of nodes visited
    GrB_Index nv = 0 ;
    LAGRAPH_OK (GrB_Vector_nvals (&nv, v)) ;
    fprintf (stderr, "number of nodes visited: %g out of %g"
        " (%g %% of the graph)\n", (double) nv, (double) n,
        100. * (double) nv / (double) n) ;

    // TODO: you can't typecast from vectors to matrices ...  (except in
    // SuiteSparse:GraphBLAS).  Need a function LAGraph_Vector_isequal.
    bool isequal = false ;
    LAGRAPH_OK (LAGraph_isequal (&isequal, (GrB_Matrix) v, (GrB_Matrix) v3,
        NULL)) ;

    if (!isequal)
    {
        fprintf (stderr, "ERROR! simple and push-pull differ\n") ;
        abort ( ) ;
    }

    LAGRAPH_OK (LAGraph_isequal (&isequal, (GrB_Matrix) v, (GrB_Matrix) v2,
        NULL)) ;
    if (!isequal)
    {
        fprintf (stderr, "ERROR! simple and method2   differ\n") ;
        abort ( ) ;
    }

    LAGRAPH_OK (LAGraph_isequal (&isequal, (GrB_Matrix) v, (GrB_Matrix) v4,
        NULL)) ;
    if (!isequal)
    {
        fprintf (stderr, "ERROR! simple and PULL   differ\n") ;
        abort ( ) ;
    }

    LAGRAPH_OK (LAGraph_isequal (&isequal, (GrB_Matrix) v, (GrB_Matrix) v5,
        NULL)) ;
    if (!isequal)
    {
        fprintf (stderr, "ERROR! simple and best   differ\n") ;
        abort ( ) ;
    }

    #if 0
    // diff = v2 - v
    LAGRAPH_OK (GrB_Vector_new (&diff, GrB_INT32, n)) ;
    LAGRAPH_OK (GrB_eWiseAdd (diff, NULL, NULL, GrB_MINUS_INT32, v2, v, NULL)) ;
    // err = or (diff)
    bool err ;
    LAGRAPH_OK (GrB_reduce (&err, NULL, LAGraph_LOR_MONOID, diff, NULL)) ;
    if (err)
    {
        GxB_fprint (diff, 3, stderr) ;
        GxB_fprint (v2, 3, stderr) ;
        GxB_fprint (v,  3, stderr) ;
        fprintf (stderr, "ERROR! simple and method2 differ\n") ;
        abort ( ) ;
    }
    else
    {
        fprintf (stderr, "results of simple and method2 are the same\n") ;
    }
    #endif

    //--------------------------------------------------------------------------
    // write the result to stdout (check them outside of this main program)
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; i < n ; i++)
    {
        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        int64_t x = 0 ;
        LAGRAPH_OK (GrB_Vector_extractElement (&x, v, i)) ;
        printf ("%" PRIu64 "\n", x) ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_finalize ( )) ;
    fprintf (stderr, "bfs_test: all tests passed\n") ;
    return (GrB_SUCCESS) ;
}

