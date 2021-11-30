//------------------------------------------------------------------------------
// LAGraph_lcc: local clustering coefficient
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

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
// structure of A is used.  This step takes extra time and memory to sanitize the
// input matrix A.  For a fair comparison in the LDBC benchmark, sanitize
// should be false.

// Results are undefined if sanitize is false, and the matrix A has any entries
// not equal to 1 (even zero-weight edges are not allowed), or if it has self
// edges.

#define LAGraph_FREE_ALL            \
{                                   \
    GrB_free (&C) ;                 \
    GrB_free (&CL) ;                \
    if (sanitize) GrB_free (&S) ;   \
    GrB_free (&U) ;                 \
    GrB_free (&W) ;                 \
    GrB_free (&LCC) ;               \
    GrB_free (&LAGraph_COMB_DIR_FP64) ;                 \
    GrB_free (&LAGraph_COMB_UNDIR_FP64) ;                 \
}

#include <LAGraph.h>
#include <LAGraphX.h>

//------------------------------------------------------------------------------

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

// z = x * (x - 1), used by LAGraph_lcc.
// This operator calculates the 2-permutation of d(v).
void LAGraph_comb_dir_fp64
(
    void *z,
    const void *x
)
{
    double xd = *(double *) x ;
    double *zd = (double *) z ;
    (*zd) = ((xd) * (xd - 1)) ;
}

// z = x * (x - 1) / 2, used by LAGraph_lcc.
// This operator calculates the 2-combination of d(v).
void LAGraph_comb_undir_fp64
(
    void *z,
    const void *x
)
{
    double xd = *(double *) x ;
    double *zd = (double *) z ;
    (*zd) = ((xd) * (xd - 1)) / 2;
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_lcc            // compute lcc for all nodes in A
(
    GrB_Vector *LCC_handle,     // output vector
    GrB_Type   *LCC_type,
    const GrB_Matrix A,         // input matrix
    bool symmetric,             // if true, the matrix is symmetric
    bool sanitize,              // if true, ensure A is binary
    double t [2]                // t [0] = sanitize time, t [1] = lcc time,
                                // in seconds
)
{
#if !(LG_SUITESPARSE)
    // FIXME: use GrB_select and make this pure GrB*
    return GrB_PANIC;
#else

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (LCC_handle == NULL || LCC_type == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    GrB_Matrix C = NULL, CL = NULL, S = NULL, U = NULL ;
    GrB_Vector W = NULL, LCC = NULL ;
    (*LCC_type) = NULL;
    GrB_UnaryOp LAGraph_COMB_DIR_FP64 = NULL ;
    GrB_UnaryOp LAGraph_COMB_UNDIR_FP64 = NULL ;
    GrB_Info info ;

    // n = size of A (# of nodes in the graph)
    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
#if LG_SUITESPARSE
    GxB_Format_Value fmt ;
    LAGRAPH_OK (GxB_get (A, GxB_FORMAT, &fmt)) ;
    if (fmt != GxB_BY_ROW)
    {
        LAGRAPH_ERROR ("A must be stored by row", GrB_INVALID_VALUE) ;
    }
#endif

    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    double tic [2] ;
    t [0] = 0 ;         // sanitize time
    t [1] = 0 ;         // LCC time

    // FIXME: use operators that ignore the values of A
    if (sanitize)
    {
        LAGraph_Tic (tic, NULL) ;

        // S = binary structure of A
        GrB_Matrix_new (&S, GrB_FP64, n, n) ;
        GrB_apply (S, NULL, NULL, GrB_ONEB_FP64, A, 0, NULL) ;

        // remove all self edges
        GrB_select (S, NULL, NULL, GrB_OFFDIAG, S, 0, NULL) ;
        LAGraph_Toc (&t[0], tic, NULL) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary with no self edges.
        // Results are undefined if this condition does not hold.
        S = A ;
    }

    LAGraph_Tic (tic, NULL) ;

    //--------------------------------------------------------------------------
    // create the operators for LAGraph_lcc
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_COMB_DIR_FP64,
                                 F_UNARY (LAGraph_comb_dir_fp64),
                                 GrB_FP64, GrB_FP64)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_COMB_UNDIR_FP64,
                                 F_UNARY (LAGraph_comb_undir_fp64),
                                 GrB_FP64, GrB_FP64)) ;

    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_Matrix_new (&U, GrB_UINT32, n, n)) ;

    if (symmetric)
    {
        C = S ;
        S = NULL ;

        //----------------------------------------------------------------------
        // U = triu(C)
        //----------------------------------------------------------------------

        LAGRAPH_OK (GxB_select (U, NULL, NULL, GxB_TRIU, C, NULL, NULL)) ;

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
        // U = triu(D)
        //----------------------------------------------------------------------

        // note that L=U' since D is symmetric
        LAGRAPH_OK (GxB_select (U, NULL, NULL, GxB_TRIU, D, NULL, NULL)) ;
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

    // CL<C> = C*L = C*U' using a masked dot product
    LAGRAPH_OK (GrB_Matrix_new (&CL, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_mxm (CL, C, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C, U,
                         GrB_DESC_T1));
    GrB_free (&U) ; U = NULL;

    //--------------------------------------------------------------------------
    // Calculate LCC
    //--------------------------------------------------------------------------

    // LCC(i) = sum (CL (i,:)) = # of triangles at each node
    LAGRAPH_OK (GrB_Vector_new (&LCC, GrB_FP64, n)) ;
    LAGRAPH_OK (GrB_reduce (LCC, NULL, NULL, GrB_PLUS_FP64, CL, NULL)) ;
    GrB_free (&CL) ; CL = NULL;

    // LCC = LCC ./ W
    LAGRAPH_OK (GrB_eWiseMult (LCC, NULL, NULL, GrB_DIV_FP64, LCC, W, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*LCC_handle) = LCC ; LCC = NULL ;
    (*LCC_type) = GrB_FP64;

    LAGraph_FREE_ALL ;
    LAGraph_Toc (&t[1], tic, NULL) ;
    return (GrB_SUCCESS) ;
#endif
}
