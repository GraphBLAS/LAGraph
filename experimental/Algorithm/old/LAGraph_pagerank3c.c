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

//------------------------------------------------------------------------------
// LAGraph_pagerank3c: GAP-style PageRank, with import/export

// See also LAGraph_pagerank3a, for the same computation without import/export.

// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619 which assumes that both A and A' are
// already available, as are the row and column degrees.

// The GAP Benchmark algorithm assumes the graph has no nodes with no out-going
// edges (otherwise, a divide-by-zero occurs).  In terms of the adjacency
// matrix, it assumes there are no rows in A that have no entries.

// For fastest results, the input matrix should stored in GxB_BY_COL format.
// TODO: or use AT by row, since the GAP assumes both A and A' are available.

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#define LAGRAPH_FREE_WORK           \
{                                   \
    LAGRAPH_FREE (I) ;              \
    LAGRAPH_FREE (pr) ;             \
    LAGRAPH_FREE (prior) ;          \
    GrB_free (&v) ;                 \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    GrB_free (result) ;             \
}

GrB_Info LAGraph_pagerank3c // PageRank definition
(
    GrB_Vector *result,     // output: array of LAGraph_PageRank structs
    GrB_Matrix A,           // binary input graph, not modified
    const float *LA_RESTRICT d_out, // out degree of each node (GrB_FP32, size n)
    float damping,          // damping factor (typically 0.85)
    int itermax,            // maximum number of iterations
    int *iters              // output: number of iterations taken
)
{

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Index n, ncols ;
    GrB_Vector v = NULL ;
    GrB_Index *I = NULL ;
    float *LA_RESTRICT pr = NULL ;
    float *prior = NULL ;
    (*result) = NULL ;

    LAGr_Matrix_ncols (&ncols, A) ;
    LAGr_Matrix_nrows (&n, A) ;

    if (ncols != n)
    {
        LAGRAPH_ERROR ("matrix must be square", GrB_DIMENSION_MISMATCH) ;
    }

    // Teleport value
    const float teleport = (1 - damping) / n ;
    const float tol = 1e-4 ;
    float rdiff = 1 ;       // first iteration is always done

    GrB_Type type = GrB_FP32 ;

    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (n, nthreads) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    // initializing pr and I
    pr = LAGraph_malloc (n, sizeof (float)) ;
    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    // do not need I
    #else
    I =  LAGraph_malloc (n, sizeof (GrB_Index)) ;
    #endif
    prior = LAGraph_malloc (n, sizeof (float)) ;
    if (pr == NULL || I == NULL || prior == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }

    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int64_t k = 0 ; k < n ; k++)
    {
        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
        // do not need I
        #else
        I [k] = k ;
        #endif
        pr [k] = 1.0/n ;
    }

    //--------------------------------------------------------------------------
    // pagerank iterations
    //--------------------------------------------------------------------------

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {
// printf ("\n============================ pagerank 3C iter: %d\n", (*iters)) ;
        // Importance calculation
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t i = 0 ; i < n; i++)
        {
            prior [i] = pr [i] ;
            pr [i] = damping * pr [i] / d_out [i] ;
        }

        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
        LAGr_Vector_import_Full (&v, GrB_FP32, n, (void **) (&pr), NULL) ;
        #else
        // import pr and I into v
        LAGr_Vector_import (&v, GrB_FP32, n, n, &I, (void **) (&pr), NULL) ;
        #endif

        // Calculate total PR of all inbound vertices: v = A' * v
        LAGr_mxv (v, NULL, NULL, GxB_PLUS_SECOND_FP32, A, v, LAGraph_desc_tooo);

        GrB_Index nvals ;
        LAGr_Vector_nvals (&nvals, v) ;
        if (nvals != n)
        {
            LAGRAPH_ERROR ("Matrix must not have empty rows or columns!",
                GrB_PANIC) ;
        }

        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
        LAGr_Vector_export_Full (&v, &type, &n, (void **) (&pr), NULL) ;
        #else
        // export v to pr and I
        LAGr_Vector_export (&v, &type, &n, &nvals, &I, (void **) (&pr), NULL) ;
        #endif

        // add teleport and check for convergence
        rdiff = 0 ;
        #pragma omp parallel for num_threads(nthreads) schedule(static) \
            reduction(+:rdiff)
        for (int64_t i = 0 ; i < n; i++)
        {
            pr [i] += teleport ;
            rdiff += fabsf (prior [i] - pr [i]) ;
        }
    }

    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    LAGr_Vector_import_Full (result, GrB_FP32, n, (void **) (&pr), NULL) ;
    #else
    // import result (pr and I) into final result
    LAGr_Vector_import (result, GrB_FP32, n, n, &I, (void **) (&pr), NULL) ;
    #endif
    LAGRAPH_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

