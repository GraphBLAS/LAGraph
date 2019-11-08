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
}

int main (int argc, char **argv)
{
    GrB_Info info;

    GrB_Matrix A_in = NULL;
    GrB_Matrix A = NULL;
    GrB_Vector path_lengths = NULL;

    LAGRAPH_OK (LAGraph_init());

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // read in the file in Matrix Market format
    LAGRAPH_OK (LAGraph_mmread(&A_in, stdin));

    // convert to boolean, pattern-only
    LAGRAPH_OK (LAGraph_pattern(&A, A_in, GrB_FP64));

    // finish any pending computations
    GrB_Index nvals;
    GrB_Matrix_nvals (&nvals, A);

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols;
    LAGr_Matrix_nrows(&nrows, A);
    LAGr_Matrix_ncols(&ncols, A);
    GrB_Index n = nrows;

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    fprintf (stderr, "\n=========="
        "input graph: nodes: %llu edges: %llu\n", n, nvals) ;

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

    LAGr_Vector_new(&path_lengths, GrB_FP64, n);

    for (int trial = 0; trial < ntrials; trial++)
    {
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
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());

    return (GrB_SUCCESS);
}

