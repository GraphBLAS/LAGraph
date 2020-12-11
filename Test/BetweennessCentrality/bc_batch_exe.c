//------------------------------------------------------------------------------
// bc_batch_exe: read in matrix and compute betweenness centrality batch version
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
// bc_batch_exe < in_file > out_file
// in_file is the Matrix Market file of the adjacency matrix.
// out_file is the betweenness centrality of all vertices.

#include "bc_test.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A);                  \
    GrB_free (&Abool);              \
    GrB_free (&v_batch);            \
}

int main (int argc, char **argv)
{
    GrB_Info info;

    GrB_Matrix A = NULL;
    GrB_Matrix Abool = NULL;
    GrB_Vector v_batch = NULL;

    LAGRAPH_OK (LAGraph_init());

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
    LAGr_Matrix_nrows(&nrows, A);
    LAGr_Matrix_ncols(&ncols, A);
    GrB_Index n = nrows;

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    printf ( "\n=========="
        "input graph: nodes: %"PRIu64" edges: %"PRIu64"\n", n, nvals);

    int nthreads = LAGraph_get_nthreads();
    printf ( "Starting bc_batch_exe\n");
    printf ( " - nthreads: %d\n", nthreads);

    //--------------------------------------------------------------------------
    // Compute betweenness centrality using batch algorithm from all nodes
    //--------------------------------------------------------------------------

    printf ( " - Start: Betweenness Centrality (Batch Algorithm)\n");

    // Create batch of vertices to use in traversal
    // int n_batch = /* size_of_batch */;
    // GrB_Index *vertex_list = malloc(sizeof(GrB_Index) * n);

    // ... or use GrB_ALL for all vertices
    int n_batch = n;
    const GrB_Index *vertex_list = GrB_ALL;

    GrB_free (&v_batch);
    LAGRAPH_OK (LAGraph_bc_batch (&v_batch, A, vertex_list, n_batch));

    printf ( " - End: Betweenness Centrality (Batch Algorithm)\n");

    //--------------------------------------------------------------------------
    // write the result to stdout
    //--------------------------------------------------------------------------

    for (int64_t i = 0; i < n; i++)
    {
        // if the entry v(i) is not present, x is unmodified, so '0' is printed
        float x = 0;
        LAGr_Vector_extractElement (&x, v_batch, i);
        printf("%f\n", x);
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());

    return (GrB_SUCCESS);
}

