//------------------------------------------------------------------------------
// sssp_test: read in a matrix and test delta stepping sssp
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

// Contributed by Jinhao Chen, Scott Kolodziej and Tim Davis, Texas A&M
// University

// usage:
// sssp_test < in.mtx > out
// in.mtx is the Matrix Market file, out is the level set.

// sssp_test in.mtx source.mtx delta > out
// sssp_test in.grb source.mtx delta > out

#include "sssp_test.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A_in);               \
    GrB_free (&A);                  \
    GrB_free (&SourceNodes) ;       \
    LAGRAPH_FREE (I);               \
    LAGRAPH_FREE (J);               \
    LAGRAPH_FREE (W);               \
    LAGRAPH_FREE (d);               \
    LAGRAPH_FREE (pi);              \
    GrB_free (&path_lengths);       \
    GrB_free (&path_lengths1);      \
}

int main (int argc, char **argv)
{
    GrB_Info info;

    GrB_Index s = 0;
    int32_t delta = 3;

    GrB_Matrix A_in = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix SourceNodes = NULL;
    GrB_Index *I = NULL, *J = NULL; // for col/row indices of entries from A
    int32_t *W = NULL;
    int32_t *d = NULL;              // for BF result
    int64_t *pi = NULL;
    GrB_Vector path_lengths = NULL;
    GrB_Vector path_lengths1 = NULL;

    bool test_pass = true;
    double tic[2];

    LAGRAPH_OK (LAGraph_init());
    //GxB_set(GxB_CHUNK, 4096) ;
    // LAGraph_set_nthreads (1) ;

    LAGraph_tic (tic);

    int batch_size = 4 ;

    //--------------------------------------------------------------------------
    // get the matrix
    //--------------------------------------------------------------------------

    if (argc > 1)
    {

        //----------------------------------------------------------------------
        // Usage: ./sssp_test matrixfile.mtx ...
        //----------------------------------------------------------------------

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
            LAGRAPH_OK (LAGraph_binread (&A_in, filename)) ;
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
            LAGRAPH_OK (LAGraph_mmread(&A_in, f));
            fclose (f) ;
        }
    }
    else
    {

        //----------------------------------------------------------------------
        // Usage:  ./sssp_test < matrixfile.mtx
        //----------------------------------------------------------------------

        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_mmread(&A_in, stdin));
    }

    // get the size of the problem.
    GrB_Index nvals;
    GrB_Matrix_nvals (&nvals, A_in);
    GrB_Index nrows, ncols;
    LAGr_Matrix_nrows(&nrows, A_in);
    LAGr_Matrix_ncols(&ncols, A_in);
    GrB_Index n = nrows;

    //--------------------------------------------------------------------------
    // get delta
    //--------------------------------------------------------------------------

    if (argc > 2)
    {
        // usage:  ./sssp_test matrix delta ...
        delta = atoi (argv [2]) ;
    }
    else
    {
        // usage:  ./sssp_test matrix
        // or:     ./sssp_test < matrix
        delta = 2 ;
    }
    printf ("delta: %d\n", delta) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    if (argc > 3)
    {

        //----------------------------------------------------------------------
        // usage:  ./sssp_test matrix delta sourcenodes
        //----------------------------------------------------------------------

        // read in source nodes in Matrix Market format from the input file
        char *filename = argv [3] ;
        printf ("sources: %s\n", filename) ;
        FILE *f = fopen (filename, "r") ;
        if (f == NULL)
        {
            printf ("Source node file not found: [%s]\n", filename) ;
            exit (1) ;
        }
        LAGRAPH_OK (LAGraph_mmread (&SourceNodes, f)) ;
        fclose (f) ;

    }
    else
    {

        //----------------------------------------------------------------------
        // create random source nodes
        //----------------------------------------------------------------------

        #define NSOURCES 1

        LAGRAPH_OK (GrB_Matrix_new (&SourceNodes, GrB_INT64, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            int64_t i = rand ( ) % n ;
            // SourceNodes [k] = i 
            LAGRAPH_OK (GrB_Matrix_setElement (SourceNodes, i, k, 0)) ;
        }
        GrB_Matrix_nvals (&nvals, SourceNodes);
    }

    double t_read = LAGraph_toc (tic) ;
    printf ("read time: %g sec\n", t_read) ;

    //--------------------------------------------------------------------------

    // convert input matrix to INT32
    GrB_Type type ;
    GxB_Matrix_type (&type, A) ;
    if (type == GrB_INT32)
    {
        A = A_in ;
        A_in = NULL ;
    }
    else
    {
        GrB_Matrix_new (&A, GrB_INT32, n, n) ;
        GrB_apply (A, NULL, NULL, GrB_IDENTITY_INT32, A_in, NULL) ;
        GrB_free (&A_in) ;
    }

    // get the number of source nodes
    GrB_Index nsource;
    LAGRAPH_OK (GrB_Matrix_nrows (&nsource, SourceNodes));
    GrB_Index ignore;
    LAGr_Matrix_nvals (&ignore, SourceNodes);

    // try converting to column format (this is slower than the default)
    // GxB_set (A, GxB_FORMAT, GxB_BY_COL) ;

    // GxB_fprint (A, 2, stdout) ;
    // GxB_fprint (SourceNodes, GxB_COMPLETE, stdout) ;

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    int nthreads = LAGraph_get_nthreads();
    printf ("input graph: nodes: %"PRIu64" edges: %"PRIu64" nthreads %d\n",
        n, nvals, nthreads) ;

    int ntrials =  1;//(int) nsource;   // TODO for GAP use nsource
    double total_time1 = 0 ;
    double total_time2 = 0 ;
    double total_time3 = 0 ;
    double total_time31 = 0 ;
    double total_time32 = 0 ;
    double t1, t2, t3;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {

        //----------------------------------------------------------------------
        // get the source node for this trial
        //----------------------------------------------------------------------

        s = -1 ;
        LAGRAPH_OK (GrB_Matrix_extractElement (&s, SourceNodes, trial, 0));
        // convert from 1-based to 0-based
        s-- ;
        //printf ("\nTrial %d : source node: %"PRIu64"\n", trial, s) ;

        //----------------------------------------------------------------------
        // Compute shortest path using delta stepping with given node and delta
        //----------------------------------------------------------------------

        #if 0
        printf(" - Start Test: delta-stepping Single Source Shortest Paths"
            " (apply operator)\n");

        // start the timer
        LAGraph_tic (tic) ;

        GrB_free (&path_lengths);
        LAGRAPH_OK (LAGraph_sssp (&path_lengths, A, s, delta)) ;

        // stop the timer
        t2 = LAGraph_toc (tic) ;
        printf ("SSSP (apply)    time: %12.6g (sec), rate:"
            " %12.6g (1e6 edges/sec)\n", t2, 1e-6*((double) nvals) / t2) ;

        total_time2 += t2;
        #endif

        //----------------------------------------------------------------------
        // Compute shortest path using delta stepping with given node and delta
        //----------------------------------------------------------------------
#if 0
        // Start the timer
        LAGraph_tic (tic);

        GrB_free (&path_lengths1);
        LAGRAPH_OK (LAGraph_sssp1 (&path_lengths1, A, s, delta)) ;

        // Stop the timer
        t3 = LAGraph_toc (tic);
        printf ("SSSP1 (select)  time: %12.6g (sec), rate:"
            " %12.6g (1e6 edges/sec)\n", t3, 1e-6*((double) nvals) / t3) ;

        total_time3 += t3;
#endif

#if 0
        // Start the timer
        LAGraph_tic (tic);
        GrB_free (&path_lengths1);
        LAGRAPH_OK (LAGraph_sssp11a(&path_lengths1, A, s, delta, true)) ;

        // Stop the timer
        t3 = LAGraph_toc (tic);
        printf ("SSSP11a(select)  time: %12.6g (sec), rate:"
            " %12.6g (1e6 edges/sec)\n", t3, 1e-6*((double) nvals) / t3) ;

        // Start the timer
        LAGraph_tic (tic);

        GrB_free (&path_lengths1);
        LAGRAPH_OK (LAGraph_sssp11b (&path_lengths1, A, s, delta, true)) ;

        // Stop the timer
        t3 = LAGraph_toc (tic);
        printf ("SSSP11b (select)  time: %12.6g (sec), rate:"
            " %12.6g (1e6 edges/sec)\n", t3, 1e-6*((double) nvals) / t3) ;
#endif


        // Start the timer
        LAGraph_tic (tic);
        GrB_free (&path_lengths1);
        LAGRAPH_OK (LAGraph_sssp11 (&path_lengths1, A, s, delta, true)) ;

        // Stop the timer
        t3 = LAGraph_toc (tic);
        //printf ("SSSP11 (select)  time: %12.6g (sec), rate:"
        //    " %12.6g (1e6 edges/sec)\n", t3, 1e-6*((double) nvals) / t3) ;

        total_time31 += t3;

#if 0
        // Start the timer
        LAGraph_tic (tic);

        GrB_free (&path_lengths1);
        LAGRAPH_OK (LAGraph_sssp2 (&path_lengths1, A, s, delta)) ;

        // Stop the timer
        t3 = LAGraph_toc (tic);
        //printf ("SSSP2 (select)  time: %12.6g (sec), rate:"
        //    " %12.6g (1e6 edges/sec)\n", t3, 1e-6*((double) nvals) / t3) ;
        total_time32 += t3;
#endif

        //----------------------------------------------------------------------
        // find shortest path using BF on node s with LAGraph_pure_c
        //----------------------------------------------------------------------

#if 0
        // get the triplet form for the Bellman-Ford function
        I = LAGraph_malloc (nvals, sizeof(GrB_Index)) ;
        J = LAGraph_malloc (nvals, sizeof(GrB_Index)) ;
        W = LAGraph_malloc (nvals, sizeof(int32_t)) ;
        if (I == NULL || J == NULL || W == NULL)
        {
            LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
        }
        LAGRAPH_OK (GrB_Matrix_extractTuples_INT32(I, J, W, &nvals, A));

        // start the timer
        LAGraph_tic (tic) ;

        LAGRAPH_FREE (d) ;
        LAGRAPH_FREE (pi) ;
        LAGRAPH_OK (LAGraph_BF_pure_c (&d, &pi, s, n, nvals, I, J, W)) ;

        // stop the timer
        t1 = LAGraph_toc (tic) ;
        printf ("BF_pure_c       time: %12.6g (sec), rate:"
            " %g (1e6 edges/sec)\n", t1, 1e-6*((double) nvals) / t1) ;

        total_time1 += t1;

        LAGRAPH_FREE (pi) ;
        LAGRAPH_FREE (I) ;
        LAGRAPH_FREE (J) ;
        LAGRAPH_FREE (W) ;
#endif

        //----------------------------------------------------------------------
        // write the result to result file if there is none
        //----------------------------------------------------------------------

#if 0
        //if( access( fname, F_OK ) == -1 )// check if the result file exists
        {
            FILE *file = fopen(fname)
            for (int64_t i = 0; i < n; i++)
            {
                fprintf(file, "%d\n", d[i]);
            }
        }
#endif

        //----------------------------------------------------------------------
        // check the result for correctness
        //----------------------------------------------------------------------

#if 0
        for (int64_t i = 0; i < n; i++)
        {
            bool test_result ;

            #if 0
            double x = INT32_MAX;
            LAGr_Vector_extractElement(&x, path_lengths, i);
            test_result = ((double)d[i] == x);
            test_pass &= test_result;
            if (!test_result)
            {
                printf ("  Failure at index %"PRId64" calculated by sssp\n", i);
                printf ("  x = %g\n", x);
                printf ("  d = %d\n", d[i]);
                printf ("\n") ;
            }
            #endif

            int32_t x1 = INT32_MAX;
            LAGr_Vector_extractElement(&x1, path_lengths1, i);
            test_result = (d[i] == x1);
            test_pass &= test_result;
            if (!test_result)
            {
                printf ("  Failure at index %"PRId64" caculated by sssp1\n", i);
                printf ("  x = %d\n", x1);
                printf ("  d = %d\n", d[i]);
                printf ("\n") ;
            }
        }
#endif

        LAGRAPH_FREE (d) ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    //printf ("ntrials: %d\n", ntrials) ;
    #if 0
    printf ("Average time per trial (Bellman-Ford pure C): %12.6g sec\n",
        total_time1 / ntrials);
    #endif
    #if 0
    printf ("Average time per trial (apply operator): %g sec\n",
        total_time2 / ntrials);
    #endif
    printf ("Average time per trial (SSSP1, with select):  %12.6g sec "
        "(delta %d)\n", total_time3 / ntrials, delta);
    printf ("Average time per trial (SSSP11, with select):  %12.6g sec "
        "(delta %d)\n", total_time31 / ntrials, delta);
    printf ("Average time per trial (SSSP2, with select):  %12.6g sec "
        "(delta %d)\n", total_time32 / ntrials, delta);

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());
#if 0
    if(!test_pass)
    {
        printf("ERROR! TEST FAILURE\n") ;
    }
    else
    {
        printf("all tests passed\n");
    }
#endif
    return (GrB_SUCCESS);
}

