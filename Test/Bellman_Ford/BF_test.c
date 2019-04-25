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
    GrB_free (&A) ;                 \
    GrB_free (&d) ;                 \
    GrB_free (&pi) ;                \
    GrB_free (&h) ;                 \
    GrB_free (&d1) ;                \
}

int main (int argc, char **argv)
{
    GrB_Info info;
    GrB_Matrix A = NULL ;
    GrB_Vector d = NULL ;
    GrB_Vector pi = NULL ;
    GrB_Vector h = NULL ;
    GrB_Vector d1 = NULL ;

    LAGRAPH_OK (LAGraph_init ( )) ;
    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // read in the file in Matrix Market format
    LAGRAPH_OK (LAGraph_mmread (&A, stdin)) ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nvals ;
    GrB_Matrix_nvals (&nvals, A) ;
    GrB_Index nrows, ncols ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;

    //--------------------------------------------------------------------------
    // set the diagonal to 0
    //--------------------------------------------------------------------------
    for (GrB_Index i = 0; i < n; i++)
    {
        LAGRAPH_OK (GrB_Matrix_setElement_FP64 (A, 0, i, i));
    }

    //--------------------------------------------------------------------------
    // get the source node
    //--------------------------------------------------------------------------

    GrB_Index s = 0 ;
    if (argc > 1)
    {
        sscanf (argv [1], "%" PRIu64, &s) ; 
    }

    fprintf (stderr, "\n=========="
        "input graph: nodes: %lu edges: %lu source node: %lu\n", n, nvals, s) ;

    //--------------------------------------------------------------------------
    // run the LAGraph_BF_full on node s
    //--------------------------------------------------------------------------

    int ntrials = 1 ;       // increase this to 10, 100, whatever, for more
                            // accurate timing

    // start the timer
    double tic [2] ;
    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d) ;
        LAGRAPH_OK (LAGraph_BF_full (&d, &pi, &h, A, s)) ;
    }

    // stop the timer
    double t1 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "FB_full   time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t1, 1e-6*((double) nvals) / t1) ;

    //--------------------------------------------------------------------------
    // run the BFS on node s with LAGraph_BF_basic
    //--------------------------------------------------------------------------

    // start the timer
    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d1) ;
        LAGRAPH_OK (LAGraph_BF_basic (&d1, A, s)) ;
    }

    // stop the timer
    double t2 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "FB_basic  time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t2, 1e-6*((double) nvals) / t2) ;
    fprintf (stderr, "speedup of FB_basic:   %g\n", t1/t2) ;


    //--------------------------------------------------------------------------
    // check results
    //--------------------------------------------------------------------------
    bool isequal = false, ok = true ; 

    LAGRAPH_OK (LAGraph_Vector_isequal (&isequal, d1, d, NULL)) ;
    if (!isequal)
    {
        fprintf (stderr, "ERROR! BF_full and BF_basic   differ\n") ;
        ok = false ;
    }

    //--------------------------------------------------------------------------
    // write the result to stdout (check them outside of this main program)
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; i < n ; i++)
    {
        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        int64_t x = 0 ;
        LAGRAPH_OK (GrB_Vector_extractElement (&x, d, i)) ;
        printf ("%" PRIu64 "\n", x) ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_finalize ( )) ;
    fprintf (stderr, "bfs_test: ") ;
    if (ok)
    {
        fprintf (stderr, "all tests passed\n") ;
    }
    else
    {
        fprintf (stderr, "TEST FAILURE\n") ;
    }
    fprintf (stderr,
    "------------------------------------------------------------\n\n") ;
    return (GrB_SUCCESS) ;
}
