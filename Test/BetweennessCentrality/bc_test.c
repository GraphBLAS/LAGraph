//------------------------------------------------------------------------------
// bc_test: read in a matrix and test betweenness centrality
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

// Contributed by Scott Kolodziej, Texas A&M University

// usage:
// bc_test < in > out
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

    bool tests_pass = true;

    LAGRAPH_OK (LAGraph_init ());

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // read in the file in Matrix Market format
    LAGRAPH_OK (LAGraph_mmread(&A, stdin));

    // convert to boolean, pattern-only
    LAGRAPH_OK (LAGraph_pattern(&Abool, A, GrB_BOOL));

    GrB_free (&A);
    A = Abool;
    Abool = NULL ;

    // finish any pending computations
    GrB_Index nvals;
    GrB_Matrix_nvals (&nvals, A);

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols;
    LAGRAPH_OK (GrB_Matrix_nrows(&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols(&ncols, A));
    GrB_Index n = nrows;

    // AT = A'
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

    GrB_Matrix_nvals (&nvals, AT);

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    fprintf (stderr, "\n=========="
        "input graph: nodes: %"PRIu64" edges: %"PRIu64"\n", n, nvals) ;

    int nthreads = LAGraph_get_nthreads();
    printf("Starting Betweenness Centrality Tests\n");
    printf(" - nthreads: %d\n", nthreads);


    int ntrials = 1;        // increase this to 10, 100, whatever, for more
                            // accurate timing

    printf(" - ntrials: %d\n", ntrials);

    //--------------------------------------------------------------------------
    // Compute betweenness centrality from all nodes (Brandes)
    //--------------------------------------------------------------------------

    printf(" - Start Test: Betweenness Centrality (Brandes Algorithm)\n");

    // Start the timer
    double tic [2];
    LAGraph_tic (tic);

    LAGRAPH_OK (GrB_Vector_new(&v_brandes, GrB_FP32, n));

    for (int trial = 0; trial < ntrials; trial++)
    {
        for (int vertex = 0; vertex < n; vertex++)
        {
            GrB_free (&v) ;
            LAGRAPH_OK (LAGraph_bc (&v, A, vertex)) ;
            LAGRAPH_OK (GrB_eWiseAdd(v_brandes, GrB_NULL, GrB_NULL, GrB_PLUS_FP32, v_brandes, v, GrB_NULL));
        }
    }

    // Stop the timer
    double t1 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "Brandes  time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t1, 1e-6*((double) nvals) / t1) ;

    printf(" - End Test: Betweenness Centrality (Brandes Algorithm)\n");

    //--------------------------------------------------------------------------
    // Compute betweenness centrality using batch algorithm from all nodes
    //--------------------------------------------------------------------------

    printf(" - Start Test: Betweenness Centrality (Batch Algorithm)\n");

    // Create batch of vertices to use in traversal
    // int n_batch = /* size_of_batch */;
    // GrB_Index *vertex_list = malloc(sizeof(GrB_Index) * n);

    // ... or use GrB_ALL for all vertices
    int n_batch = n;
    const GrB_Index *vertex_list = GrB_ALL;

    // Start the timer
    LAGraph_tic (tic) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&v_batch) ;
//      LAGRAPH_OK (LAGraph_bc_batch (&v_batch, A, vertex_list, n_batch)) ;
        LAGRAPH_OK (LAGraphX_bc_batch3 (&v_batch, A, AT, vertex_list, n_batch));
    }

    // Stop the timer
    double t2 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "Batch    time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t2, 1e-6*((double) nvals) / t2) ;

    printf(" - End Test: Betweenness Centrality (Batch Algorithm)\n");

    //--------------------------------------------------------------------------
    // write the result to stdout
    //--------------------------------------------------------------------------

    printf(" - Betweenness Centrality Numerical Results\n\n");

    printf("   +-------------------------+\n");
    printf("   | v_i | Brandes |  Batch  |\n");
    printf("   +-------------------------+\n");

    for (int64_t i = 0; i < n; i++)
    {
        printf("   | %3"PRId64" ", i);

        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        float x1 = 0;
        LAGRAPH_OK (GrB_Vector_extractElement (&x1, v_brandes, i));
        printf("| %7.2f ", x1);

        float x2 = 0;
        LAGRAPH_OK (GrB_Vector_extractElement (&x2, v_batch, i));
        printf ("| %7.2f |\n", x2);

        // Check that both methods give the same results
        bool test_result = (fabs(x1 - x2) / (1E-10 + fmax(x1, x2)) < 1E-5);
        tests_pass &= test_result;
        if (!test_result)
        {
            fprintf(stderr, "Failure at index %"PRId64"\n", i);
            fprintf(stderr, "x1 = %f\n", x1);
            fprintf(stderr, "x2 = %f\n", x2);
            fprintf(stderr, "Error = %f\n", fabs(x1-x2) / (1E-6 + fmax(x1,x2)));
        }
    }

    printf("   +-------------------------+\n");

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());
    fprintf (stderr, "bc_test: ");
    if (tests_pass)
    {
        fprintf (stderr, "all tests passed\n");
        printf("all tests passed\n");
    }
    else
    {
        fprintf (stderr, "TEST FAILURE\n");
        printf("TEST FAILURE\n");
    }
    fprintf (stderr,
    "------------------------------------------------------------\n\n");
    return (GrB_SUCCESS);
}

