//------------------------------------------------------------------------------
// p3test: read in (or create) a matrix and test PageRank3
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

// Contributed by Tim Davis, Texas A&M and Gabor Szarnyas, BME

// usage:
// p3test < in > out

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL                        \
{                                               \
    if (P != NULL) { free (P) ; P = NULL ; }    \
    GrB_free (&A) ;                             \
    GrB_free (&PR) ;                            \
    GrB_free (&A_temp) ;                        \
}

int main (int argc, char **argv)
{

    GrB_Info info ;
    GrB_Matrix A = NULL ;
    GrB_Matrix A_temp = NULL ;
    LAGraph_PageRank *P = NULL ;
    GrB_Vector PR = NULL;

    LAGRAPH_OK (LAGraph_init ( )) ;

    int nthreads_max = LAGraph_get_nthreads ( ) ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;

    if (argc > 1)
    {
        // Usage:
        //      ./p3test matrixfile.mtx
        //      ./p3test matrixfile.grb

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
            FILE *f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_mmread(&A, f));
            fclose (f) ;
        }

    }
    else
    {

        // Usage:  ./p3test < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_mmread(&A, stdin));
    }


    // GxB_fprint (A, GxB_SHORT, stdout) ;
    // LAGraph_mmwrite (A, stdout) ;

    // convert to FP32
    LAGRAPH_OK (LAGraph_pattern (&A_temp, A, GrB_BOOL)) ;
    // LAGraph_mmwrite (A_temp, stdout) ;
    GrB_free (&A) ;
    A = A_temp ;
    A_temp = NULL ;
    LAGRAPH_OK(GxB_set (A, GxB_FORMAT, GxB_BY_COL));
    // GxB_fprint (A, GxB_COMPLETE, stdout) ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;

    //--------------------------------------------------------------------------
    // ensure the matrix has no empty rows
    //--------------------------------------------------------------------------

    GrB_Vector d = NULL ;
    LAGRAPH_OK(GrB_Vector_new(&d, GrB_FP32, n));
    LAGRAPH_OK(GrB_reduce( d, NULL, NULL, GxB_PLUS_FP32_MONOID, A, NULL ));
    GrB_Index non_dangling ;
    LAGRAPH_OK (GrB_Vector_nvals (&non_dangling, d)) ;

    if (non_dangling < n)
    {
        // A = A+I if A has any dangling nodes
        printf ("Matrix has %"PRId64" empty rows\n", n - non_dangling) ;
        for (GrB_Index i = 0; i < n; i++)
        {
            float di = 0 ;
            LAGRAPH_OK (GrB_Vector_extractElement (&di, d, i)) ;
            if (di == 0)
            {
                non_dangling++ ;
                LAGRAPH_OK (GrB_Matrix_setElement (A, 1, i, i));
            }
        }
    }

    if (non_dangling < n) { printf ("ERROR!\n") ; abort ( ) ; }
    GrB_free (&d) ;

    // finish any pending computations
    GrB_Index nvals ;
    GrB_Matrix_nvals (&nvals, A) ;

    // LAGRAPH_OK (GrB_Matrix_setElement (A, 0, 0, n-1)) ;     // hack

    printf ("\n=========="
            "input graph: nodes: %"PRIu64" edges: %"PRIu64"\n", n, nvals) ;

    double tread = LAGraph_toc (tic) ;
    printf ("read time: %g sec\n", tread) ;

    GxB_fprint (A, GxB_SHORT, stdout) ;

    //--------------------------------------------------------------------------
    // compute the pagerank (both methods)
    //--------------------------------------------------------------------------

    int ntrials = 16 ;       // increase this to 10, 100, whatever, for more
    // accurate timing

    float tol = 1e-4 ;
    int iters, itermax = 100 ;

    // #define NTHRLIST 7
    // int nthread_list [NTHRLIST] = {1, 2, 4, 8, 16, 20, 40} ;

    // #define NTHRLIST 2
    // int nthread_list [NTHRLIST] = {1, 40} ;    

    #define NTHRLIST 1
    int nthread_list [NTHRLIST] = {40} ;    

    nthread_list [NTHRLIST-1] = nthreads_max ;

    //--------------------------------------------------------------------------
    // method 3a
    //--------------------------------------------------------------------------

    for (int kk = 0 ; kk < NTHRLIST; kk++)
    {
        int nthreads = nthread_list [kk] ;
        LAGraph_set_nthreads (nthreads) ;
        
        double total_time = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            GrB_free (&PR) ;
            LAGraph_tic (tic) ;
            LAGRAPH_OK (LAGraph_pagerank3a (&PR, A, 0.85, itermax, &iters)) ;
            double t1 = LAGraph_toc (tic) ;
            total_time += t1 ;
            printf ("trial %2d, time %16.6g\n", trial, t1) ;
        }

        double t = total_time / ntrials ;
        printf ("Average pagerank3a  time: %16.6g (sec), "
                "rate: %10.4g (1e6 edges/sec) iters: %d threads: %d\n",
                t, 1e-6*((double) nvals) / t, iters, nthreads) ;

    }

    GxB_Vector_fprint(PR, "---- PR ------", GxB_SHORT, stdout);
    GrB_free (&PR) ;

    //--------------------------------------------------------------------------
    // method 3c
    //--------------------------------------------------------------------------

    for (int kk = 0 ; kk < NTHRLIST; kk++)
    {
        int nthreads = nthread_list [kk] ;
        LAGraph_set_nthreads (nthreads) ;

        double total_time = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            GrB_free (&PR) ;
            LAGraph_tic (tic) ;
            LAGRAPH_OK (LAGraph_pagerank3c (&PR, A, 0.85, itermax, &iters)) ;
            double t1 = LAGraph_toc (tic) ;
            total_time += t1 ;
            printf ("trial %2d, time %16.6g\n", trial, t1) ;
        }

        double t = total_time / ntrials ;
        printf ("Average pagerank3c  time: %16.6g (sec), "
                "rate: %10.4g (1e6 edges/sec) iters: %d threads: %d\n",
                t, 1e-6*((double) nvals) / t, iters, nthreads) ;
    }

    GxB_Vector_fprint(PR, "---- PR ------", GxB_SHORT, stdout);

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_finalize ( )) ;
    return (GrB_SUCCESS) ;
}
