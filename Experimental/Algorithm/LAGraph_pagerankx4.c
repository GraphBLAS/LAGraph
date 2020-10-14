//------------------------------------------------------------------------------
// LAGraph_pagerankx4: pagerank using a real semiring
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2020 LAGraph Contributors.

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
// LAGraph_pagerankx4: GAP-style PageRank, with import/export

// Tim Davis and Mohsen Aznaveh.

// See also LAGraph_pagerank3f, for the same computation without import/export.
// This version is just slightly faster than LAGraph_pagerank3f (perhaps 10%
// at most, sometimes the difference is smaller).

// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619 which assumes that both A and A' are
// already available, as are the row and column degrees.

// The GAP Benchmark algorithm assumes the graph has no nodes with no out-going
// edges (otherwise, a divide-by-zero occurs when dividing by d_out [i] below).
// In terms of the adjacency matrix, it assumes there are no rows in A that
// have no entries.

// For fastest results, the input matrix should stored in GxB_BY_COL format.
// TODO: or use AT by row, since the GAP assumes both A and A' are available.

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)

    // no need for vi and wi
    #define LAGRAPH_FREE_WORK           \
    {                                   \
        LAGRAPH_FREE (vx) ;             \
        LAGRAPH_FREE (wx) ;             \
        LAGRAPH_FREE (prior) ;          \
        GrB_free (&v) ;                 \
        GrB_free (&w) ;                 \
    }

#else

    #define LAGRAPH_FREE_WORK           \
    {                                   \
        LAGRAPH_FREE (vi) ;             \
        LAGRAPH_FREE (vx) ;             \
        LAGRAPH_FREE (wi) ;             \
        LAGRAPH_FREE (wx) ;             \
        LAGRAPH_FREE (prior) ;          \
        GrB_free (&v) ;                 \
        GrB_free (&w) ;                 \
    }

#endif

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    GrB_free (result) ;             \
}

GrB_Info LAGraph_pagerankx4 // PageRank definition
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
    GrB_Index n ;
    GrB_Vector v = NULL, w = NULL ;
    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    // no need for vi and wi
    #else
    GrB_Index *vi = NULL, *wi = NULL ;
    #endif
    float *LA_RESTRICT vx = NULL ;
    float *LA_RESTRICT wx = NULL ;
    float *LA_RESTRICT prior = NULL ;
    GrB_Type type = GrB_FP32 ;
    (*result) = NULL ;
    LAGr_Matrix_nrows (&n, A) ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,2,0) )
    GrB_Descriptor desc_t0 = GrB_DESC_T0 ;
    #else
    GrB_Descriptor desc_t0 = LAGraph_desc_tooo ;
    #endif

    const float teleport = (1 - damping) / n ;
    const float tol = 1e-4 ;
    float rdiff = 1 ;       // first iteration is always done

    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (n, nthreads) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    // allocate workspace
    vx = LAGraph_malloc (n, sizeof (float)) ;
    wx = LAGraph_malloc (n, sizeof (float)) ;
    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    // no need for vi and wi
    #else
    vi = LAGraph_malloc (n, sizeof (GrB_Index)) ;
    wi = LAGraph_malloc (n, sizeof (GrB_Index)) ;
    if (vi == NULL || wi == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }
    #endif
    prior = LAGraph_malloc (n, sizeof (float)) ;
    if (vx == NULL || prior == NULL || wx == NULL)
    {
        LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
    }

    // v = 1/n
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int64_t k = 0 ; k < n ; k++)
    {
        vx [k] = 1.0 / n ;
        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
        // no need for vi and wi
        #else
        vi [k] = k ;
        wi [k] = k ;
        #endif
    }

    //--------------------------------------------------------------------------
    // pagerank iterations
    //--------------------------------------------------------------------------

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {
        // prior = v ;
        // v = damping * v ./ dout ;
        // w (:) = teleport
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t i = 0 ; i < n ; i++)
        {
            prior [i] = vx [i] ;
            vx [i] = damping * vx [i] / d_out [i] ;
            wx [i] = teleport ;
        }

        // import wx and wi into w
        // import vx and vi into v
        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
        LAGr_Vector_import_Full (&w, type, n, (void **) (&wx), NULL) ;
        LAGr_Vector_import_Full (&v, type, n, (void **) (&vx), NULL) ;
        #else
        LAGr_Vector_import (&w, type, n, n, &wi, (void **) (&wx), NULL) ;
        LAGr_Vector_import (&v, type, n, n, &vi, (void **) (&vx), NULL) ;
        #endif

        // w += A'*v
        LAGr_mxv (w, NULL, GrB_PLUS_FP32, GxB_PLUS_SECOND_FP32, A, v, desc_t0) ;

        // export w to vx and vi (the new score; note the swap)
        // export v to wx and wi (workspace for next iteration)
        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
        LAGr_Vector_export_Full (&w, &type, &n, (void **) (&vx), NULL) ;
        LAGr_Vector_export_Full (&v, &type, &n, (void **) (&wx), NULL) ;
        #else
        LAGr_Vector_export (&w, &type, &n, &n, &vi, (void **) (&vx), NULL) ;
        LAGr_Vector_export (&v, &type, &n, &n, &wi, (void **) (&wx), NULL) ;
        #endif

        // check for convergence
        rdiff = 0 ;
        #pragma omp parallel for num_threads(nthreads) schedule(static) \
            reduction(+:rdiff)
        for (int64_t i = 0 ; i < n ; i++)
        {
            rdiff += fabsf (prior [i] - vx [i]) ;
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    LAGr_Vector_import_Full (result, type, n, (void **) (&vx), NULL) ;
    #else
    LAGr_Vector_import (result, type, n, n, &vi, (void **) (&vx), NULL) ;
    #endif
    LAGRAPH_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

