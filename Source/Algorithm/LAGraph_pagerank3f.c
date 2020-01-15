//------------------------------------------------------------------------------
// LAGraph_pagerank3f: pagerank using a real semiring
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
// LAGraph_pagerank3f: GAP-style PageRank, all work done in GraphBLAS

// See also LAGraph_pagerank3c, for the same computation but with import/export.

// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619 which assumes that both A and A' are
// already available, as are the row and column degrees.

// The GAP Benchmark algorithm assumes the graph has no nodes with no out-going
// edges (otherwise, a divide-by-zero occurs).  In terms of the adjacency
// matrix, it assumes there are no rows in A that have no entries.

// For fastest results, the input matrix A should be stored in GxB_BY_COL
// format (TODO: check this on input).  The values of A are ignored; just its
// pattern is used.  All nodes of A must have at least one out-going edge
// (that is, the matrix A cannot have any empty rows); otherwise, a
// divide-by-zero occurs and the results are undefined.

// Contributed by Tim Davis and Mohsen Aznaveh.

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGr_free (&r) ;                \
    LAGr_free (&d) ;                \
    LAGr_free (&t) ;                \
    LAGr_free (&w) ;                \
}

GrB_Info LAGraph_pagerank3f // PageRank definition
(
    GrB_Vector *result,     // output: array of LAGraph_PageRank structs
    GrB_Matrix A,           // binary input graph, not modified
    GrB_Vector d_out,       // outbound degree of all nodes (not modified)
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
    GrB_Vector r = NULL, d = NULL, t = NULL, w = NULL ;

    LAGr_Matrix_nrows (&n, A) ;

    // r = 1 / n
    LAGr_Vector_new (&t, GrB_FP32, n) ;
    LAGr_Vector_new (&r, GrB_FP32, n) ;
    LAGr_Vector_new (&w, GrB_FP32, n) ;
    LAGr_assign (r, NULL, NULL, 1.0 / n, GrB_ALL, n, NULL) ;

    const float teleport = (1 - damping) / n ;

    // prescale with damping factor, so it isn't done each iteration
    // d = d_out / damping ;
    LAGr_Vector_dup (&d, d_out) ;
    LAGr_assign (d, NULL, GrB_DIV_FP32, damping, GrB_ALL, n, NULL) ;

    const float tol = 1e-4 ;
    float rdiff = 1 ;       // so first iteration is always done

    //--------------------------------------------------------------------------
    // pagerank iterations
    //--------------------------------------------------------------------------

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {
        // swap t and r ; now t is the old score
        GrB_Vector temp = t ; t = r ; r = temp ;

        // w = t ./ d
        LAGr_eWiseMult (w, NULL, NULL, GrB_DIV_FP32, t, d, NULL) ;

        // r = teleport
        LAGr_assign (r, NULL, NULL, teleport, GrB_ALL, n, NULL) ;

        // r += A'*w
        LAGr_mxv (r, NULL, GrB_PLUS_FP32, GxB_PLUS_SECOND_FP32, A, w,
            LAGraph_desc_tooo);

        // t -= r
        LAGr_assign (t, NULL, GrB_MINUS_FP32, r, GrB_ALL, n, NULL) ;

        // t = abs (t)
        LAGr_apply (t, NULL, NULL, GxB_ABS_FP32, t, NULL) ;

        // rdiff = sum (t)
        LAGr_reduce (&rdiff, NULL, GxB_PLUS_FP32_MONOID, t, NULL) ;
   }

    (*result) = r ;
    r = NULL ;
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

