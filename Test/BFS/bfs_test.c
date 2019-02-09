//------------------------------------------------------------------------------
// bfs_test: read in (or create) a matrix and test BFS
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// usage:
// bfs_test s < in > out
// s is the staring node, in is the Matrix Market file, out is the level set.

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&AT) ;                \
    GrB_free (&A) ;                 \
    GrB_free (&Abool) ;             \
    GrB_free (&v) ;                 \
    GrB_free (&v2) ;                \
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

    //--------------------------------------------------------------------------
    // get the source node
    //--------------------------------------------------------------------------

    GrB_Index s = 0 ;
    if (argc > 1)
    {
        sscanf (argv [1], "%" PRIu64, &s) ; 
    }

    //--------------------------------------------------------------------------
    // run the BFS on node s
    //--------------------------------------------------------------------------

    // start the timer
    double tic [2] ;
    LAGraph_tic (tic) ;

    LAGRAPH_OK (LAGraph_bfs_simple (&v, A, s, n)) ;

    // stop the timer
    double t1 = LAGraph_toc (tic) ;
    fprintf (stderr, "simple    time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t1, 1e-6*((double) nvals) / t1) ;

    //--------------------------------------------------------------------------
    // now the BFS on node s using push-pull instead
    //--------------------------------------------------------------------------

    // AT = A'
    LAGRAPH_OK (GrB_Matrix_new (&AT, GrB_BOOL, ncols, nrows)) ;
    LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;

    LAGraph_tic (tic) ;

    LAGRAPH_OK (LAGraph_bfs_pushpull (&v2, A, AT, s, n)) ;

    double t2 = LAGraph_toc (tic) ;
    fprintf (stderr, "push/pull time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t2, 1e-6*((double) nvals) / t2) ;
    fprintf (stderr, "speedup of push/pull: %g\n", t1/t2) ;

    //--------------------------------------------------------------------------
    // check results
    //--------------------------------------------------------------------------

    // TODO: you can't typecast from vectors to matrices ...  (except in
    // SuiteSparse:GraphBLAS).  Need a function LAGraph_Vector_isequal.
    bool isequal = false ;
    LAGRAPH_OK (LAGraph_isequal (&isequal, (GrB_Matrix) v, (GrB_Matrix) v2,
        NULL)) ;

    if (isequal)
    {
        fprintf (stderr, "results of simple and push-pull are the same\n") ;
    }
    else
    {
        fprintf (stderr, "ERROR! simple and push-pull differ\n") ;
        abort ( ) ;
    }

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

