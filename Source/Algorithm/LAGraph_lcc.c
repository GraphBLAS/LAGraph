//------------------------------------------------------------------------------
// LAGraph_lcc.c: local clustering coefficient
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_lcc.c authors:  Gabor Szarnyas and Balint Hegyi, TU Budapest
// (with accented characters: G\'{a}bor Sz\'{a}rnyas and B\'{a}lint Hegyi,
// using LaTeX syntax). https://inf.mit.bme.hu/en/members/szarnyasg

// This function was originally written for the LDBC Graphalytics benchmark,
// at https://graphalytics.org/ .

// Modified and added to LAGraph by Tim Davis, with permission by the authors.

// The local clustering coefficient is a measure for each node of a directed
// graph.  It is fully described in the following document:
// https://ldbc.github.io/ldbc_graphalytics_docs/graphalytics_spec.pdf

// For each node v, the lcc(v) is the ratio between the number of edges between
// neighbors of the node v, and the maximum possible number of edges between
// these neighbors.  If a node has fewer than 2 neighbors, then its coefficient
// is defined as zero (TODO where is this done?).

// Let N_in(v)  = the set of nodes u such that (u,v) is an edge.
// Let N_out(v) = the set of nodes u such that (v,u) is an edge.
// Let N(v) = union (N_in(v), N_out(v)).
// Then the metric lcc(v) is defined as:

// lcc(v) = (sum for all u in N(v) of |intersection (N(v), N_out(u))) /
//          ( |N(v)| * (|N(v)|-1))

// That is, for directed graphs, the set of neighbors N(v) is found without
// taking directions into account, but a node u that has both an edge (u,v) and
// (v,u) is counted just once.  However, edges directions are enforced when
// considering two nodes u1 and u2 that are both in N(v).

// The input matrix A must be square.  If A is known to be binary (with all
// explicit edge weights equal to 1), then sanitize can be false.  This is
// the case for the LDBC benchmark.

// Otherwise, if sanitize is true, edge weights of A are ignored and only the
// pattern of A is used.  This step takes extra time and memory to sanitize the
// input matrix A.  For a fair comparison in the LDBC benchmark, sanitize
// should be false.

// Results are undefined if A has non-binary edge and sanitize is false.

// TODO what about self-edges?  They should be ignored, I assume?

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&LCC) ;               \
    GrB_free (&AT) ;                \
    if (sanitize) GrB_free (&S) ;   \
    GrB_free (&C) ;                 \
    GrB_free (&CA) ;                \
    GrB_free (&W) ;                 \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_lcc            // compute lcc for all nodes in A
(
    GrB_Vector *LCC_handle,     // output vector
    const GrB_Matrix A,         // input matrix
    bool sanitize               // if true, ensure A is binary
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (LCC_handle == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    GrB_Matrix AT = NULL, C = NULL, CA = NULL, S = NULL ;
    GrB_Vector W = NULL, LCC = NULL ;
    GrB_Info info ;

    // n = size of A (# of nodes in the graph)
    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    //--------------------------------------------------------------------------
    // ensure input is binary
    //--------------------------------------------------------------------------

    if (sanitize)
    {
        // S = binary pattern of A
        LAGRAPH_OK (LAGraph_pattern (&S, A)) ;
        // TODO remove self edges here?
    }
    else
    {
        // use the input as-is, assume it is binary
        S = A ;
    }

    //--------------------------------------------------------------------------
    // C = A+A' to create an undirected graph C
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Matrix_new (&AT, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, S, NULL)) ;
    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_eWiseAdd (C, NULL, NULL, GrB_LOR, S, AT, NULL)) ;
    if (sanitize) GrB_free (&S) ;
    GxB_print (A, 3) ;
    GxB_print (C, 3) ;

    //--------------------------------------------------------------------------
    // Find wedges of each node
    //--------------------------------------------------------------------------

    // W(i) = sum (C (i,:))
    LAGRAPH_OK (GrB_Vector_new (&W, GrB_FP64, n)) ;
    LAGRAPH_OK (GrB_reduce (W, NULL, NULL, GrB_PLUS_FP64, C, NULL)) ;

    // Create vector W for containing number of wedges per vertex
    // W(i) = W(i) * (W(i)-1)
    LAGRAPH_OK (GrB_apply (W, NULL, NULL, LAGraph_COMB_FP64, W, NULL)) ;
    GxB_print (W, 3) ;

    //--------------------------------------------------------------------------
    // Calculate triangles
    //--------------------------------------------------------------------------

    // CA<C> = C*A = C*AT' using a masked dot product
    LAGRAPH_OK (GrB_Matrix_new (&CA, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_mxm (CA, C, NULL, LAGraph_PLUS_TIMES_FP64, C, AT,
        LAGraph_desc_otoo)) ;
    GrB_free (&AT) ;

    // CA = CA.*C via element-wise multiplication
    LAGRAPH_OK (GrB_eWiseMult (CA, NULL, NULL, GrB_TIMES_FP64, CA, C, NULL)) ;
    GrB_free (&C) ;

    //--------------------------------------------------------------------------
    // Calculate LCC
    //--------------------------------------------------------------------------

    // LCC(i) = sum (CA (i,:)) = # of triangles at each node
    LAGRAPH_OK (GrB_Vector_new (&LCC, GrB_FP64, n)) ;
    LAGRAPH_OK (GrB_reduce (LCC, NULL, NULL, GrB_PLUS_FP64, CA, NULL)) ;
    GrB_free (&CA) ;

    // LCC = LCC ./ W
    LAGRAPH_OK (GrB_eWiseMult (LCC, NULL, NULL, GrB_DIV_FP64, LCC, W, NULL)) ;
    GxB_print (LCC, 3) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*LCC_handle) = LCC ;
    LCC = NULL ;            // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

