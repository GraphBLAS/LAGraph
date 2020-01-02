//------------------------------------------------------------------------------
// bc_gap_test: betweenness centrality for the GAP benchmark
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

// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University

// usage:
// bc_gap_test matrixfile.mtx sourcenodes.mtx
// in is the Matrix Market file, out is the level set.

#include "bc_test.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A);                  \
    GrB_free (&AT);                 \
    GrB_free (&Abool);              \
    GrB_free (&v);                  \
    GrB_free (&v_brandes);          \
    GrB_free (&v_batch);            \
    GrB_free (&SourceNodes) ;       \
}

int main (int argc, char **argv)
{
    GrB_Info info;
    uint64_t seed = 1;

    GrB_Matrix A = NULL;
    GrB_Matrix AT = NULL;
    GrB_Matrix Abool = NULL;
    GrB_Vector v = NULL;
    GrB_Vector v_brandes = NULL;
    GrB_Vector v_batch = NULL;
    GrB_Matrix SourceNodes = NULL ;

    bool tests_pass = true;

    LAGRAPH_OK (LAGraph_init ());
    // GxB_set (GxB_NTHREADS, 1) ;
    // GxB_set (GxB_CHUNK, 1) ;

    // Start the timer
    double tic [2];
    LAGraph_tic (tic);

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    int batch_size = 4 ;

    if (argc > 1)
    {
        // Usage:
        //      ./bc_gap_test matrixfile.mtx sourcenodes.mtx
        //      ./bc_gap_test matrixfile.grb sourcenodes.mtx

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

        // read in source nodes in Matrix Market format from the input file
        filename = argv [2] ;
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

        // Usage:  ./bc_gap_test < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_mmread(&A, stdin));

        // use nodes [1 2 3 4...] as the source nodes (in 1-based notation)
        LAGRAPH_OK (GrB_Matrix_new (&SourceNodes, GrB_INT64, batch_size, 1)) ;
        for (int64_t i = 0 ; i < batch_size ; i++)
        {
            LAGRAPH_OK (GrB_Matrix_setElement (SourceNodes, i+1, i, 0)) ;
        }
    }

    double t_read = LAGraph_toc (tic) ;
    printf ("read time: %g sec\n", t_read) ;

    LAGraph_tic (tic);

    // convert to pattern-only
    LAGRAPH_OK (LAGraph_pattern(&Abool, A, GrB_BOOL));

    GrB_free (&A);
    A = Abool;
    Abool = NULL ;

    // finish any pending computations
    GrB_Index nvals;
    GrB_Matrix_nvals (&nvals, SourceNodes);
    GrB_Matrix_nvals (&nvals, A);

    GxB_fprint (A, 2, stdout) ;

    GxB_fprint (SourceNodes, GxB_COMPLETE, stdout) ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols;
    LAGRAPH_OK (GrB_Matrix_nrows(&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols(&ncols, A));
    GrB_Index n = nrows;

    GrB_Index nsource ;
    LAGRAPH_OK (GrB_Matrix_nrows(&nsource, SourceNodes));
    if (nsource % batch_size != 0)
    {
        printf ("SourceNode size must be multiple of batch_size (%d)\n",
            batch_size) ;
        exit (1) ;
    }

    double t_setup = LAGraph_toc (tic) ;
    printf ("setup time: %g sec\n", t_setup) ;

    // AT = A'
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
        GxB_fprint (AT, 2, stdout) ;
    }
    double t_transpose = LAGraph_toc (tic) ;
    printf ("transpose time: %g\n", t_transpose) ;

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    printf ("\n========== input graph: nodes: %"PRIu64" edges: %"PRIu64"\n",
        n, nvals) ;

    int nthreads = LAGraph_get_nthreads();
    printf("TESTING LAGraphX_bc_batch3 (saxpy in both phases, nthreads %d\n",
        nthreads) ;

    int ntrials = 0 ;
    double total_time_1 = 0 ;
    double total_time_2 = 0 ;

    for (int64_t kstart = 0 ; kstart <  nsource ; kstart += batch_size)
    {

        //----------------------------------------------------------------------
        // Create batch of vertices to use in traversal
        //----------------------------------------------------------------------

        ntrials++ ;
        printf ("\nTrial %d : sources: [", ntrials) ;
        GrB_Index vertex_list [batch_size] ;
        for (int64_t k = 0 ; k < batch_size ; k++)
        {
            // get the kth source node
            GrB_Index source = -1 ;
            LAGRAPH_OK (GrB_Matrix_extractElement (&source, SourceNodes, k + kstart, 0)) ;
            // subtract one to convert from 1-based to 0-based
            source-- ;
            vertex_list [k] = source  ;
            printf (" %"PRIu64, source) ;
        }
        printf (" ]\n") ;

        //----------------------------------------------------------------------
        // Compute betweenness centrality from four source nodes (Brandes)
        //----------------------------------------------------------------------

        #if 0
        // Start the timer
        LAGraph_tic (tic);

        LAGRAPH_OK (GrB_Vector_new(&v_brandes, GrB_FP64, n));
        for (int64_t k = 0 ; k < batch_size ; k++)
        {
            // get the kth source node
            GrB_Index source = vertex_list [k] ;
            LAGRAPH_OK (LAGraph_bc (&v, A, source)) ;
//          LAGRAPH_OK (LAGraph_bc2 (&v, A, source)) ;
            LAGRAPH_OK (GrB_eWiseAdd(v_brandes, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, v_brandes, v, GrB_NULL));
            GrB_free (&v) ;
        }

        // Stop the timer
        double t1 = LAGraph_toc (tic) ;
        printf ("Brandes  time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
            t1, 1e-6*((double) nvals) / t1) ;

        total_time_1 += t1 ;
        #endif

        //----------------------------------------------------------------------
        // Compute betweenness centrality using batch algorithm
        //----------------------------------------------------------------------

        // Start the timer
        LAGraph_tic (tic) ;

//      LAGRAPH_OK (LAGraph_bc_batch  (&v_batch, A, vertex_list, batch_size)) ;
//      LAGRAPH_OK (LAGraphX_bc_batch (&v_batch, A, vertex_list, batch_size)) ;
//      LAGRAPH_OK (LAGraphX_bc_batch2 (&v_batch, A, vertex_list, batch_size)) ;
        LAGRAPH_OK (LAGraphX_bc_batch3 (&v_batch, A, AT, vertex_list, batch_size)) ;

#if 0
        LAGRAPH_OK (GrB_Vector_new(&v_batch, GrB_FP64, n));
        for (int64_t k = 0 ; k < 4 ; k++)
        {
            // get the kth source node
            GrB_Index source = vertex_list [k] ;
//          LAGRAPH_OK (LAGraph_bc (&v, A, source)) ;
            LAGRAPH_OK (LAGraph_bc2 (&v, A, source)) ;
            LAGRAPH_OK (GrB_eWiseAdd(v_batch, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, v_batch, v, GrB_NULL));
            GrB_free (&v) ;
        }
#endif

        // Stop the timer
        double t2 = LAGraph_toc (tic) ;
        printf ("Batch    time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
            t2, 1e-6*((double) nvals) / t2) ;

        total_time_2 += t2 ;

        //----------------------------------------------------------------------
        // check result
        //----------------------------------------------------------------------

        #if 0
        for (int64_t i = 0; i < n; i++)
        {

            // if the ith entry is not present, x is unmodified, so '0'
            // is printed
            double x1 = 0;
            LAGRAPH_OK (GrB_Vector_extractElement (&x1, v_brandes, i));

            double x2 = 0;
            LAGRAPH_OK (GrB_Vector_extractElement (&x2, v_batch, i));

            double err = (fabs(x1 - x2) / (1E-10 + fmax(x1, x2))) ;

            // Check that both methods give the same results
            bool test_result =  err < 1E-5;
            tests_pass &= test_result;
            if (!test_result)
            {
                printf ("%"PRId64 " %g %g %g ", i, x1, x2, err) ;
                printf ("FAIL") ;
                /*
                printf ("   Failure at index %"PRId64"\n", i);
                printf ("   x1 = %g\n", x1);
                printf ("   x2 = %g\n", x2);
                printf ("   Error = %f\n", fabs(x1-x2) / (1E-6 + fmax(x1,x2)));
                */
                // break ;
                printf ("\n") ;
            }
        }
        #endif

        GrB_free (&v_brandes) ;
        GrB_free (&v_batch) ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    printf ("ntrials: %d\n", ntrials) ;
    if (total_time_1 > 0)
    {
        printf ("Average time per trial (Brandes): %g sec\n",
            total_time_1 / ntrials);
    }
    printf ("Average time per trial (batch):   %g sec\n",
        total_time_2 / ntrials);

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());
    if (tests_pass)
    {
        printf("%s: all tests passed\n", argv [0]);
    }
    else
    {
        printf("%s: TEST FAILURE\n", argv [0]);
    }
    return (GrB_SUCCESS);
}

