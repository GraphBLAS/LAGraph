//------------------------------------------------------------------------------
// bc_exe: read in a matrix and compute betweenness centrality
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
// bc_exe < in_file > out_file
// in_file is the Matrix Market file of the adjacency matrix.
// out_file is the betweenness centrality of all vertices.

#include "bc_test.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A);                  \
    GrB_free (&Abool);              \
    GrB_free (&v);                  \
    GrB_free (&v_brandes);          \
}

int main (int argc, char **argv)
{
    GrB_Info info;
    uint64_t seed = 1;

    GrB_Matrix A = NULL;
    GrB_Matrix Abool = NULL;
    GrB_Vector v = NULL;
    GrB_Vector v_brandes = NULL;

    LAGRAPH_OK (LAGraph_init ());

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // read in the file in Matrix Market format
    LAGRAPH_OK (LAGraph_mmread(&A, stdin));

    // convert to boolean, pattern-only
    LAGRAPH_OK (LAGraph_pattern(&Abool, A, GrB_INT64));

    GrB_free (&A);
    A = Abool;

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

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    fprintf (stderr, "\n=========="
        "input graph: nodes: %"PRIu64" edges: %"PRIu64"\n", n, nvals) ;

    int nthreads = LAGraph_get_nthreads();
    fprintf(stderr, "Starting bc_exe\n");
    fprintf(stderr, " - nthreads: %d\n", nthreads);

    //--------------------------------------------------------------------------
    // Compute betweenness centrality from all nodes (Brandes)
    //--------------------------------------------------------------------------

    fprintf(stderr, " - Start: Betweenness Centrality (Brandes Algorithm)\n");

    LAGRAPH_OK (GrB_Vector_new(&v_brandes, GrB_FP64, n));

    for (int vertex = 0; vertex < n; vertex++)
    {
        GrB_free (&v) ;
        LAGRAPH_OK (LAGraph_bc (&v, A, vertex)) ;
        LAGRAPH_OK (GrB_eWiseAdd(v_brandes, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, v_brandes, v, GrB_NULL));
    }

    fprintf(stderr, " - End: Betweenness Centrality (Brandes Algorithm)\n");

    //--------------------------------------------------------------------------
    // write the result to stdout
    //--------------------------------------------------------------------------

    for (int64_t i = 0; i < n; i++)
    {
        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        float x = 0;
        LAGRAPH_OK (GrB_Vector_extractElement (&x, v_brandes, i));
        printf("%f\n", x);
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_finalize());

    return (GrB_SUCCESS) ;
}

