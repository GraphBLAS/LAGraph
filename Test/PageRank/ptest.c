//------------------------------------------------------------------------------
// ptest: read in (or create) a matrix and test PageRank
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// usage:
// ptest < in > out

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL                        \
{                                               \
    if (P != NULL) { free (P) ; P = NULL ; }    \
    GrB_free (&A) ;                             \
    GrB_free (&Abool) ;                         \
}

int main ( )
{

    GrB_Info info ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    LAGraph_PageRank *P = NULL ;

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
    Abool = NULL ;
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

    fprintf (stderr, "\n=========="
        "input graph: nodes: %"PRIu64" edges: %"PRIu64"\n", n, nvals) ;

    //--------------------------------------------------------------------------
    // compute the pagerank
    //--------------------------------------------------------------------------

    int nthreads_max = 1 ;
    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    GxB_get (GxB_NTHREADS, &nthreads_max) ;
    #endif

    int ntrials = 1 ;       // increase this to 10, 100, whatever, for more
                            // accurate timing
    
    double tol = 1e-5 ;
    int iters, itermax = 100 ;

    for (int nthreads = 1 ; nthreads <= nthreads_max ; nthreads++)
    {
        #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
        GxB_set (GxB_NTHREADS, nthreads) ;
        #endif
        printf ("\nptest nthreads %d ======================================\n",
            nthreads) ;

        // start the timer
        double tic [2] ;
        LAGraph_tic (tic) ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            if (P != NULL) { free (P) ; P = NULL ; }
            LAGRAPH_OK (LAGraph_pagerank (&P, A, itermax, tol, &iters)) ;
        }

        // stop the timer
        double t1 = LAGraph_toc (tic) / ntrials ;
        fprintf (stderr, "pagerank  time: %12.6e (sec), "
            "rate: %g (1e6 edges/sec) iters: %d threads: %d\n",
            t1, 1e-6*((double) nvals) / t1, iters, nthreads) ;
    }

    //--------------------------------------------------------------------------
    // print results
    //--------------------------------------------------------------------------

    for (int64_t k = 0 ; k < n ; k++)
    {
        printf ("%" PRIu64 " %g\n", P [k].page, P [k].pagerank) ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_finalize ( )) ;
    return (GrB_SUCCESS) ;
}

