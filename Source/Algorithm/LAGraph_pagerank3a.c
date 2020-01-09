//------------------------------------------------------------------------------
// LAGraph_pagerank3a: pagerank using a real semiring
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
// LAGraph_pagerank3a: GAP-style PageRank, all work done in GraphBLAS

// See also LAGraph_pagerank3c, for the same computation but with import/export.

// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619 which assumes that both A and A' are
// already available, as are the row and column degrees.

// The GAP Benchmark algorithm assumes the graph has no nodes with no out-going
// edges (otherwise, a divide-by-zero occurs).  In terms of the adjacency
// matrix, it assumes there are no rows in A that have no entries.

// For fastest results, the input matrix A should be stored in GxB_BY_COL
// format (TODO: check this on input).  All entries in A must be equal to 1
// (TODO: relax this condition).  All nodes of A must have at least one
// out-going edge (that is, the matrix A cannot have any empty rows).  For
// fastest results, the matrix A should not have any empty columns.

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGr_free (&v) ;                \
    LAGr_free (&pr) ;               \
    LAGr_free (&prior) ;            \
    LAGr_free (&op_diff) ;          \
};

void ddiff (void *z, const void *x, const void *y)
{
    // z = fabs (x-y)
    float delta = (* ((float *) x)) - (* ((float *) y)) ;
    (*((float *) z)) = fabs (delta) ;
}

GrB_Info LAGraph_pagerank3a // PageRank definition
(
    GrB_Vector *result,     // output: array of LAGraph_PageRank structs
    GrB_Matrix A,           // binary input graph, not modified
    GrB_Vector d_out,       // outbound degree of all nodes
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
    GrB_Vector v = NULL ;
    GrB_Vector pr = NULL ;
    GrB_BinaryOp op_diff = NULL ;
    GrB_Vector prior = NULL ;

    LAGr_Matrix_nrows (&n, A) ;

    // pr = 1 / n
    LAGr_Vector_new (&pr, GrB_FP32, n) ;
    LAGr_assign (pr, NULL, NULL, 1.0 / n, GrB_ALL, n, NULL) ;

    LAGr_Vector_new (&v, GrB_FP32, n) ;

    // Teleport value
    const float teleport = (1 - damping) / n ;

    // create binary operator to compute z = fabs (x-y)
    LAGr_BinaryOp_new (&op_diff, ddiff, GrB_FP32, GrB_FP32, GrB_FP32) ;

    const float tol = 1e-4 ;
    float rdiff = 1 ;       // so first iteration is always done

    //--------------------------------------------------------------------------
    // pagerank iterations
    //--------------------------------------------------------------------------

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {
        // prior = pr ; deep copy
        GrB_Vector_dup (&prior, pr) ;

        // Divide prior PageRank with # of outbound edges: v = pr ./ d_out
        LAGr_eWiseMult (v, NULL, NULL, GrB_DIV_FP32, pr, d_out, NULL) ;

        // Multiply importance by damping factor: v *= damping
        LAGr_assign (v, NULL, GrB_TIMES_FP32, damping, GrB_ALL, n, NULL) ;

        // Calculate total PR of all inbound vertices: v = A' * v
        LAGr_mxv (v, NULL, NULL, GxB_PLUS_SECOND_FP32, A, v, LAGraph_desc_tooo);

        // PageRank summarization: pr = (1-df)/n
        LAGr_assign (pr, NULL, NULL, teleport, GrB_ALL, n, NULL) ;

        // pr += v
        LAGr_eWiseAdd (pr, NULL, NULL, GrB_PLUS_FP32, pr, v, NULL) ;

        // rdiff = sum (|pr-prior|)
        LAGr_eWiseAdd (prior, NULL, NULL, op_diff, prior, pr, NULL) ;
        LAGr_reduce (&rdiff, NULL, GxB_PLUS_FP32_MONOID, prior, NULL) ;
        LAGr_free (&prior) ;
   }

    (*result) = pr ;
    pr = NULL ;
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

