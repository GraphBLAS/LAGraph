//------------------------------------------------------------------------------
// LAGraph_pagerank3c: pagerank using a real semiring
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

// LAGraph_pagerank3c: Alternative PageRank implementation using a real
// semiring.
//
// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619

// This algorithm assumes the graph has no nodes with no out-going edges.  In
// terms of the adjacency matrix, it assumes there are no rows in A that have
// no entries.

// For fastest results, the input matrix should stored in GxB_BY_COL format.

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE (I) ;              \
    LAGRAPH_FREE (pr) ;             \
    LAGRAPH_FREE (oldpr) ;          \
    GrB_free (&importance_vec) ;    \
    GrB_free (&grb_pr) ;            \
};

GrB_Info LAGraph_pagerank3c // PageRank definition
(
    GrB_Vector *result,    // output: array of LAGraph_PageRank structs
    GrB_Matrix A,          // binary input graph, not modified
    const float *restrict d_out, // out degree of each node (GrB_FP32, size n)
    float damping_factor,  // damping factor
    unsigned long itermax, // maximum number of iterations
    int* iters             // output: number of iterations taken
)
{

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Index n, ncols ;
    GrB_Vector importance_vec = NULL ;
    GrB_Vector grb_pr = NULL;
    GrB_Index *I = NULL ;
    float *restrict pr = NULL ;
    float *oldpr = NULL ;

    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    GrB_Index nvals ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, A)) ;

    if (ncols != n)
    {
        LAGRAPH_ERROR ("matrix must be square", GrB_DIMENSION_MISMATCH) ;
    }

    // Teleport value
    const float teleport = (1 - damping_factor) / n;
    float tol = 1e-4 ;
    float rdiff = 1 ;       // first iteration is always done

    GrB_Type type = GrB_FP32 ;
    GrB_Index d_nvals;
    GrB_Index d_n;

    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (n, nthreads) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    //--------------------------------------------------------------------------
    // compute the pagerank
    //--------------------------------------------------------------------------

    // initializing pr and I
    pr = LAGraph_malloc (n, sizeof (float)) ;
    I =  LAGraph_malloc (n, sizeof (GrB_Index)) ;
    oldpr = LAGraph_malloc (n, sizeof (float)) ;

    if (pr == NULL || I == NULL || oldpr == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }

    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int64_t k = 0 ; k < n ; k++)
    {
        I [k] = k ;
        pr [k] = 1.0/n ;
    }

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {

        // Importance calculation
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t i = 0 ; i < n; i++)
        {
            oldpr [i] = pr [i] ;
            pr [i] = damping_factor * pr [i] / d_out [i] ;
        }

        // import pr and I into importance_vec
        LAGRAPH_OK (GxB_Vector_import (&importance_vec, GrB_FP32, n, n, &I,
            (void **) (&pr), NULL)) ;

        // Calculate total PR of all inbound vertices
        // importance_vec = A' * importance_vec
        LAGRAPH_OK (GrB_mxv (importance_vec, NULL, NULL, GxB_PLUS_SECOND_FP32,
            A, importance_vec, LAGraph_desc_tooo)) ;

        GrB_Vector_nvals (&nvals, importance_vec) ;
        if (nvals != n)
        {
            LAGRAPH_ERROR ("Matrix must not have empty rows or columns!",
                GrB_PANIC) ;
        }

        // export importance_vec to pr and I
        GrB_Index nvals_exp ;
        LAGRAPH_OK (GxB_Vector_export (&importance_vec, &type, &n, &nvals_exp,
            &I, (void **) (&pr), NULL)) ;

        // add teleport and check for convergence
        rdiff = 0 ;
        #pragma omp parallel for num_threads(nthreads) schedule(static) \
            reduction(+:rdiff)
        for (int64_t i = 0 ; i < n; i++)
        {
            pr [i] += teleport ;
            rdiff += fabsf (oldpr [i] - pr [i]) ;
        }
    }

    // import result (pr and I) into grb_pr
    LAGRAPH_OK (GxB_Vector_import (&grb_pr, GrB_FP32, n, n, &I,
        (void **) (&pr), NULL)) ;
    (*result) = grb_pr ;
    grb_pr = NULL ;

    LAGRAPH_FREE (oldpr) ;
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

