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

#include "sssp_test.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A_in);               \
    GrB_free (&A);                  \
    GrB_free (&path_lengths);       \
    GrB_free (&path_lengths1);      \
}

int main (int argc, char **argv)
{
    GrB_Info info;

    GrB_Matrix A_in = NULL;
    GrB_Matrix A = NULL;
    GrB_Vector path_lengths = NULL;
    GrB_Vector path_lengths1 = NULL;

    LAGRAPH_OK (LAGraph_init());

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // TODO copy from bc_gap_test (binary vs mtx market)

    // read in the file in Matrix Market format
    LAGRAPH_OK (LAGraph_mmread(&A_in, stdin));

    // finish any pending computations
    GrB_Index nvals;
    GrB_Matrix_nvals (&nvals, A_in);

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols;
    LAGr_Matrix_nrows(&nrows, A);
    LAGr_Matrix_ncols(&ncols, A);
    GrB_Index n = nrows;

    // convert input matrix to INT32
    GrB_Matrix_new (A, GrB_INT32, n, n) ;
    GrB_apply (A, NULL, NULL, GrB_IDENTITY_INT32, A_in, NULL) ;
    GxB_print (A, 2) ;

    GrB_free (&A_in) ;

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    fprintf (stderr, "\n=========="
        "input graph: nodes: %lu edges: %lu\n", n, nvals) ;

    int nthreads = LAGraph_get_nthreads();
    fprintf(stderr, "Starting Single Source Shortest Paths Tests\n");
    fprintf(stderr, " - nthreads: %d\n", nthreads);


    int ntrials = 1;        // increase this to 10, 100, whatever, for more
                            // accurate timing

    fprintf(stderr, " - ntrials: %d\n", ntrials);

    //--------------------------------------------------------------------------
    // Compute betweenness centrality from all nodes (Brandes)
    //--------------------------------------------------------------------------

    fprintf(stderr, " - Start Test: Single Source Shortest Paths\n");

    // Start the timer
    double tic [2];
    LAGraph_tic (tic);

    // TODO use 64 trials (or whatever SourceNodes is)

    for (int trial = 0; trial < ntrials; trial++)
    {
        GrB_free (&path_lengths);
        LAGRAPH_OK (LAGraph_sssp (&path_lengths, A, 0, 3)) ;
    }

    // Stop the timer
    double t1 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "SSSP+Delta Stepping  time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t1, 1e-6*((double) nvals) / t1) ;

    fprintf(stderr, " - End Test: Single Source Shortest Paths\n");

    //--------------------------------------------------------------------------
    // write the result to stdout
    //--------------------------------------------------------------------------

    for (int64_t i = 0; i < n; i++)
    {
        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        double x = 0;
        LAGr_Vector_extractElement(&x, path_lengths, i);
        printf("%f\n", x);
    }


    //--------------------------------------------------------------------------
    // Compute betweenness centrality from all nodes (Brandes)
    //--------------------------------------------------------------------------

    fprintf(stderr, " - Start Test: Single Source Shortest Paths\n");

    // Start the timer
    LAGraph_tic (tic);


    for (int trial = 0; trial < ntrials; trial++)
    {
        GrB_free (&path_lengths1);
        LAGRAPH_OK (LAGraph_sssp1 (&path_lengths1, A, 0, 3)) ;
    }

    // Stop the timer
    t1 = LAGraph_toc (tic) / ntrials ;
    fprintf (stderr, "SSSP+Delta Stepping  time: %12.6e (sec), rate: %g (1e6 edges/sec)\n",
        t1, 1e-6*((double) nvals) / t1) ;

    fprintf(stderr, " - End Test: Single Source Shortest Paths\n");

    //--------------------------------------------------------------------------
    // write the result to stdout
    //--------------------------------------------------------------------------

    bool test_pass = true;
    for (int64_t i = 0; i < n; i++)
    {
        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        double x = 0, x1 = 0;
        LAGr_Vector_extractElement(&x1, path_lengths1, i);
        LAGr_Vector_extractElement(&x, path_lengths, i);
	test_pass &= (x == x1);
    }
    if(!test_pass)
    {
	printf ("ERROR! TEST FAILURE\n") ;
    }
    else
    {
	printf ("all tests passed\n");
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());

    return (GrB_SUCCESS);
}

