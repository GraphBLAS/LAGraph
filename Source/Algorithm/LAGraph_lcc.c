//------------------------------------------------------------------------------
// LAGraph_lcc: local clustering coefficient
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

// LAGraph_lcc: Contributed by Gabor Szarnyas and Balint Hegyi,
// Budapest University of Technology and Economics
// (with accented characters: G\'{a}bor Sz\'{a}rnyas and B\'{a}lint Hegyi,
// using LaTeX syntax). https://inf.mit.bme.hu/en/members/szarnyasg .
// Modified by Tim Davis.

// This function was originally written for the LDBC Graphalytics benchmark,
// at https://graphalytics.org/ .

// The local clustering coefficient is a measure for each node of a graph.
// Its definition is fully described in the following document:
// https://ldbc.github.io/ldbc_graphalytics_docs/graphalytics_spec.pdf

// For each node v, the lcc(v) is the ratio between the number of edges between
// neighbors of the node v, and the maximum possible number of edges between
// these neighbors.  If a node v has fewer than 2 neighbors, then its
// coefficient is defined as zero, and the vth entry does not appear in the
// sparse vector LCC returned.

// Let N_in(v)  = the set of nodes u such that (u,v) is an edge.
// Let N_out(v) = the set of nodes u such that (v,u) is an edge.
// Let N(v) = union (N_in(v), N_out(v)).
// Then the metric lcc(v) is defined as:

// lcc(v) = (sum for all u in N(v) of |intersection (N(v), N_out(u))) /
//          ( |N(v)| * (|N(v)|-1) )

// That is, for directed graphs, the set of neighbors N(v) is found without
// taking directions into account, but a node u that has both an edge (u,v) and
// (v,u) is counted just once.  However, edge directions are enforced when
// considering two nodes u1 and u2 that are both in N(v), i.e. when counting
// the number of edges between neighbors, (u,v) and (v,u) are counted as two.
// To account for this, the maximum possible number of edges for vertex v is
// determined as the 2-combination of |N(v)| for undirected graphs and as the
// 2-permutation of |N(v)| for directed graphs.

// The input matrix A must be square.  If A is known to be binary (with all
// explicit edge weights equal to 1), then sanitize can be false.  This is the
// case for the LDBC benchmark.

// Otherwise, if sanitize is true, edge weights of A are ignored and only the
// pattern of A is used.  This step takes extra time and memory to sanitize the
// input matrix A.  For a fair comparison in the LDBC benchmark, sanitize
// should be false.

// Results are undefined if sanitize is false, and the matrix A has any entries
// not equal to 1 (even zero-weight edges are not allowed), or if it has self
// edges.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&C) ;                 \
    GrB_free (&CL) ;                \
    if (sanitize) GrB_free (&S) ;   \
    GrB_free (&L) ;                 \
    GrB_free (&W) ;                 \
    GrB_free (&LCC) ;               \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_lcc            // compute lcc for all nodes in A
(
    GrB_Vector *LCC_handle,     // output vector
    const GrB_Matrix A,         // input matrix
    bool symmetric,             // if true, the matrix is symmetric
    bool sanitize,              // if true, ensure A is binary
    double t [2]                // t [0] = sanitize time, t [1] = lcc time,
                                // in seconds
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (LCC_handle == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    GrB_Matrix C = NULL, CL = NULL, S = NULL, L = NULL ;
    GrB_Vector W = NULL, LCC = NULL ;
    GrB_Info info ;

    // n = size of A (# of nodes in the graph)
    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    double tic [2] ;
    t [0] = 0 ;         // sanitize time
    t [1] = 0 ;         // LCC time

    if (sanitize)
    {
        LAGraph_tic (tic) ;

        // S = binary pattern of A
        LAGRAPH_OK (LAGraph_pattern (&S, A)) ;

        // remove all self edges
        LAGRAPH_OK (LAGraph_prune_diag (S)) ;
        t [0] = LAGraph_toc (tic) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary with no self edges.
        // Results are undefined if this condition does not hold.
        S = A ;
    }

    LAGraph_tic (tic) ;

    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_Matrix_new (&L, GrB_UINT32, n, n)) ;

    if (symmetric)
    {
        C = S ;
        S = NULL ;

        //----------------------------------------------------------------------
        // L = tril(C)
        //----------------------------------------------------------------------
        LAGRAPH_OK (GxB_select (L, NULL, NULL, GxB_TRIL, C, NULL, NULL)) ;
    }
    else
    {
        GrB_Matrix AT = NULL, D = NULL ;

        LAGRAPH_OK (GrB_Matrix_new (&AT, GrB_FP64, n, n)) ;
        LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, S, NULL)) ;

        //----------------------------------------------------------------------
        // C = A \/ A' to create an undirected graph C
        //----------------------------------------------------------------------

        LAGRAPH_OK (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
        LAGRAPH_OK (GrB_eWiseAdd (C, NULL, NULL, GrB_LOR, S, AT, NULL)) ;

        //----------------------------------------------------------------------
        // D = A + A' to create an undirected multigraph D
        //----------------------------------------------------------------------

        LAGRAPH_OK (GrB_Matrix_new (&D, GrB_FP64, n, n)) ;
        LAGRAPH_OK (GrB_eWiseAdd (D, NULL, NULL, GrB_PLUS_FP64, S, AT, NULL)) ;

        GrB_free (&AT) ;
        if (sanitize) GrB_free (&S) ;

        //----------------------------------------------------------------------
        // L = tril(D)
        //----------------------------------------------------------------------

        LAGRAPH_OK (GxB_select (L, NULL, NULL, GxB_TRIL, D, NULL, NULL)) ;
        GrB_free (&D) ;
    }

    //--------------------------------------------------------------------------
    // Find wedges of each node
    //--------------------------------------------------------------------------

    // W(i) = sum (C (i,:))
    LAGRAPH_OK (GrB_Vector_new (&W, GrB_FP64, n)) ;
    LAGRAPH_OK (GrB_reduce (W, NULL, NULL, GrB_PLUS_FP64, C, NULL)) ;

    // Compute vector W defining the number of wedges per vertex
    if (symmetric)
    {
        // the graph is undirected
        LAGRAPH_OK (GrB_apply(W, NULL, NULL, LAGraph_COMB_UNDIR_FP64, W, NULL));
    }
    else
    {
        // the graph is directed
        LAGRAPH_OK (GrB_apply(W, NULL, NULL, LAGraph_COMB_DIR_FP64, W, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // Calculate triangles
    //--------------------------------------------------------------------------

    // CL<C> = C*L using a masked dot product
    LAGRAPH_OK (GrB_Matrix_new (&CL, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_mxm (CL, C, NULL, LAGraph_PLUS_TIMES_FP64, C, L, NULL)) ;
    GrB_free (&L) ;

    //--------------------------------------------------------------------------
    // Calculate LCC
    //--------------------------------------------------------------------------

    // LCC(i) = sum (CL (i,:)) = # of triangles at each node
    LAGRAPH_OK (GrB_Vector_new (&LCC, GrB_FP64, n)) ;
    LAGRAPH_OK (GrB_reduce (LCC, NULL, NULL, GrB_PLUS_FP64, CL, NULL)) ;
    GrB_free (&CL) ;

    // LCC = LCC ./ W
    LAGRAPH_OK (GrB_eWiseMult (LCC, NULL, NULL, GrB_DIV_FP64, LCC, W, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*LCC_handle) = LCC ;
    LCC = NULL ;            // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL ;
    t [1] = LAGraph_toc (tic) ;
    return (GrB_SUCCESS) ;
}

