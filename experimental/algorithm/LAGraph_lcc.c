//------------------------------------------------------------------------------
// LAGraph_lcc: local clustering coefficient
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabor Szarnyas and Balint Hegyi, Budapest University of
// Technology and Economics (with accented characters: G\'{a}bor Sz\'{a}rnyas
// and B\'{a}lint Hegyi, using LaTeX syntax).
// https://inf.mit.bme.hu/en/members/szarnyasg .
// Modified by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

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

#define LG_FREE_ALL                 \
{                                   \
    GrB_free (&CL) ;                \
    GrB_free (&T) ;                 \
    GrB_free (&AT) ;                \
    GrB_free (&D) ;                 \
    GrB_free (&G) ;                 \
    GrB_free (&U) ;                 \
    GrB_free (&W) ;                 \
    GrB_free (&LCC) ;               \
    GrB_free (&true_value) ;        \
    GrB_free (&zeros) ;             \
    GrB_free (&LAGraph_COMB_DIR_FP64) ;   \
    GrB_free (&LAGraph_COMB_UNDIR_FP64) ; \
}

#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

//------------------------------------------------------------------------------

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

#define LAGRAPH_COMB_DIR_FP64                                                 \
"void LAGraph_comb_dir_fp64                                               \n" \
"(                                                                        \n" \
"    void *z,                                                             \n" \
"    const void *x                                                        \n" \
")                                                                        \n" \
"{                                                                        \n" \
"    double xd = *(double *) x ;                                          \n" \
"    double *zd = (double *) z ;                                          \n" \
"    (*zd) = ((xd) * (xd - 1));                                           \n" \
"}                                                                        \n" \
"}"

#define LAGRAPH_COMB_UNDIR_FP64                                               \
"void LAGraph_comb_undir_fp64                                             \n" \
"(                                                                        \n" \
"    void *z,                                                             \n" \
"    const void *x                                                        \n" \
")                                                                        \n" \
"{                                                                        \n" \
"    double xd = *(double *) x ;                                          \n" \
"    double *zd = (double *) z ;                                          \n" \
"    (*zd) = ((xd) * (xd - 1)) / 2;                                       \n" \
"}                                                                        \n" \
"}"

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

int LAGraph_lcc            // compute lcc for all nodes in A
(
    GrB_Vector *LCC_handle,     // output vector
    const GrB_Matrix A,         // input matrix
    bool symmetric,             // if true, the matrix is symmetric
    bool sanitize,              // if true, ensure A is binary
    double t [2],               // t [0] = sanitize time, t [1] = lcc time,
                                // in seconds
    char *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    if (LCC_handle == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    GrB_Matrix AT = NULL, D = NULL, G = NULL ;
    GrB_Matrix C = NULL, CL = NULL, S = NULL, U = NULL, T = NULL ;
    GrB_Vector W = NULL, LCC = NULL ;
    GrB_Vector zeros = NULL;
    GrB_Scalar true_value = NULL;
    GrB_UnaryOp LAGraph_COMB_DIR_FP64 = NULL ;
    GrB_UnaryOp LAGraph_COMB_UNDIR_FP64 = NULL ;
    GrB_Info info ;

#if !LAGRAPH_SUITESPARSE
    LG_ASSERT (false, GrB_NOT_IMPLEMENTED) ;
#else

    // n = size of A (# of nodes in the graph)
    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
#if LAGRAPH_SUITESPARSE
    GxB_Format_Value fmt ;
    GRB_TRY (GxB_get (A, GxB_FORMAT, &fmt)) ;
    if (fmt != GxB_BY_ROW)
    {
        return (GrB_INVALID_VALUE) ;
    }
#endif

    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    t [0] = 0 ;         // sanitize time
    t [1] = 0 ;         // LCC time

    if (sanitize)
    {
        t [0] = LAGraph_WallClockTime ( ) ;

        // T = binary structure of A; note that T is iso-valued
        GrB_Matrix_new (&T, GrB_FP64, n, n) ;
        GrB_apply (T, NULL, NULL, GrB_ONEB_FP64, A, 0, NULL) ;

        // remove all self edges
        GrB_select (T, NULL, NULL, GrB_OFFDIAG, T, 0, NULL) ;
        t [0] = LAGraph_WallClockTime ( ) - t [0] ;
        S = T ;
    }
    else
    {
        // Use the input as-is, and assume it is binary with no self edges.
        // Results are undefined if this condition does not hold.
        S = A ;
    }

    // The matrix S must not be freed.  T is freed instead.  S is either A
    // (the input matrix which cannot be freed) or T (which is freed itself).
    // So the matrix S itself is not free.

    t [1] = LAGraph_WallClockTime ( ) ;

    //--------------------------------------------------------------------------
    // create the operators for LAGraph_lcc
    //--------------------------------------------------------------------------

    GRB_TRY (GxB_UnaryOp_new (&LAGraph_COMB_DIR_FP64,
        F_UNARY (LAGraph_comb_dir_fp64),
        GrB_FP64, GrB_FP64,
        "LAGraph_comb_dir_fp64",
        LAGRAPH_COMB_DIR_FP64
        )) ;

    GRB_TRY (GxB_UnaryOp_new (&LAGraph_COMB_UNDIR_FP64,
        F_UNARY (LAGraph_comb_undir_fp64),
        GrB_FP64, GrB_FP64,
        "LAGraph_comb_undir_fp64",
        LAGRAPH_COMB_UNDIR_FP64)) ;

    GRB_TRY (GrB_Matrix_new (&U, GrB_FP64, n, n)) ;

    if (symmetric)
    {

        //----------------------------------------------------------------------
        // use C as the undirected graph, and compute U = triu(C)
        //----------------------------------------------------------------------

        C = S ;     // C must not be freed
        GRB_TRY (GrB_select (U, NULL, NULL, GrB_TRIU, C, 0, NULL)) ;

    }
    else
    {

        // NOTE: AT is iso-valued
        GRB_TRY (GrB_Matrix_new (&AT, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_transpose (AT, NULL, GrB_ONEB_FP64, S, NULL)) ;

        //----------------------------------------------------------------------
        // C = A \/ A' to create an undirected graph G
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Scalar_new (&true_value, GrB_FP64)) ;
        GRB_TRY (GrB_Scalar_setElement (true_value, true)) ;

        GRB_TRY (GrB_Matrix_new (&G, GrB_FP64, n, n)) ;
        GRB_TRY (GxB_eWiseUnion (G, NULL, NULL, GrB_ONEB_FP64, S, true_value,
            AT, true_value, NULL)) ;
        GrB_free (&true_value) ;

        //----------------------------------------------------------------------
        // D = A + A' to create an undirected multigraph D
        //----------------------------------------------------------------------

        // NOTE: D is not iso-valued since it is a multigraph; it will contain
        // 1s and 2s.
        GRB_TRY (GrB_Matrix_new (&D, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_eWiseAdd (D, NULL, NULL, GrB_PLUS_FP64, S, AT, NULL)) ;
        GrB_free (&AT) ;

        //----------------------------------------------------------------------
        // U = triu(D)
        //----------------------------------------------------------------------

        // note that L=U' since D is symmetric
        // NOTE: U is not iso-valued since it is a multigraph; it will contain
        // 1s and 2s
        GRB_TRY (GrB_select (U, NULL, NULL, GrB_TRIU, D, 0, NULL)) ;
        GrB_free (&D) ;

        C = G ;     // C must not be freed; G is freed instead
    }

    //--------------------------------------------------------------------------
    // Find wedges of each node
    //--------------------------------------------------------------------------

    // W(i) = sum (C (i,:))
    GRB_TRY (GrB_Vector_new (&W, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&zeros, GrB_FP64, n)) ;
    GRB_TRY (GrB_assign (zeros, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_mxv (W, NULL, NULL, LAGraph_plus_one_fp64, C, zeros, NULL)) ;
    GrB_free (&zeros);

    // Compute vector W defining the number of wedges per vertex
    if (symmetric)
    {
        // the graph is undirected
        GRB_TRY (GrB_apply(W, NULL, NULL, LAGraph_COMB_UNDIR_FP64, W, NULL));
    }
    else
    {
        // the graph is directed
        GRB_TRY (GrB_apply(W, NULL, NULL, LAGraph_COMB_DIR_FP64, W, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // Calculate triangles
    //--------------------------------------------------------------------------

    // CL<C> = C*L = C*U' using a masked dot product
    GRB_TRY (GrB_Matrix_new (&CL, GrB_FP64, n, n)) ;

    GrB_Semiring sr = symmetric ? LAGraph_plus_one_fp64 : GrB_PLUS_TIMES_SEMIRING_FP64;
    GRB_TRY (GrB_mxm (CL, C, NULL, sr, C, U, GrB_DESC_ST1));
    GrB_free (&U) ;     // this sets U to NULL

    //--------------------------------------------------------------------------
    // Calculate LCC
    //--------------------------------------------------------------------------

    // LCC(i) = sum (CL (i,:)) = # of triangles at each node
    GRB_TRY (GrB_Vector_new (&LCC, GrB_FP64, n)) ;
    GRB_TRY (GrB_reduce (LCC, NULL, NULL, GrB_PLUS_FP64, CL, NULL)) ;
    GrB_free (&CL) ;        // this sets CL to NULL

    // LCC = LCC ./ W
    GRB_TRY (GrB_eWiseMult (LCC, NULL, NULL, GrB_DIV_FP64, LCC, W, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*LCC_handle) = LCC ; LCC = NULL ;

    LG_FREE_ALL ;
    t [1] = LAGraph_WallClockTime ( ) - t [1] ;

    return (GrB_SUCCESS) ;
#endif
}
