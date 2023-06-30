//------------------------------------------------------------------------------
// LAGraph_lcc: local clustering coefficient
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2023 by The LAGraph Contributors, All Rights Reserved.
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
// Modified by Pascal Costanza, Intel, Belgium

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
    GrB_free (&C) ;                 \
    GrB_free (&CL) ;                \
    GrB_free (&S) ;                 \
    GrB_free (&U) ;                 \
    GrB_free (&W) ;                 \
    GrB_free (&x) ;                 \
    GrB_free (&LCC) ;               \
    GrB_free (&LAGraph_COMB_FP64) ; \
}

#include "LG_internal.h"
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
"}"

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
"}"

//------------------------------------------------------------------------------

int LAGraph_lcc            // compute lcc for all nodes in A
(
    GrB_Vector *LCC_handle,     // output vector
    LAGraph_Graph G,            // input graph
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

    GrB_Matrix C = NULL, CL = NULL, S = NULL, U = NULL ;
    GrB_Vector W = NULL, LCC = NULL, x = NULL ;
    GrB_UnaryOp LAGraph_COMB_FP64 = NULL ;
    GrB_Info info ;

#if !LAGRAPH_SUITESPARSE
    LG_ASSERT (false, GrB_NOT_IMPLEMENTED) ;
#else

    GrB_Matrix A = G->A ;

    // n = size of A (# of nodes in the graph)
    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Matrix_new (&S, GrB_FP64, n, n));
    GRB_TRY (GrB_apply (S, GrB_NULL, GrB_NULL, GrB_ONEB_FP64, A, 0, GrB_NULL));
    if (G->nself_edges != 0) {
        GRB_TRY (GrB_select (S, GrB_NULL, GrB_NULL, GrB_OFFDIAG, S, 0, GrB_NULL));
    }

    //--------------------------------------------------------------------------
    // create the operators for LAGraph_lcc
    //--------------------------------------------------------------------------

    if (G->is_symmetric_structure == LAGraph_TRUE) {
        GRB_TRY (GxB_UnaryOp_new(&LAGraph_COMB_FP64,
                                 F_UNARY(LAGraph_comb_undir_fp64),
                                 GrB_FP64, GrB_FP64,
                                 "LAGraph_comb_undir_fp64",
                                 LAGRAPH_COMB_UNDIR_FP64
        ));
    } else {
        GRB_TRY (GxB_UnaryOp_new(&LAGraph_COMB_FP64,
                                 F_UNARY(LAGraph_comb_dir_fp64),
                                 GrB_FP64, GrB_FP64,
                                 "LAGraph_comb_dir_fp64",
                                 LAGRAPH_COMB_DIR_FP64
        ));
    }

    GRB_TRY (GrB_Matrix_new (&U, GrB_FP64, n, n)) ;

    if (G->is_symmetric_structure == LAGraph_TRUE)
    {
        C = S ;
        S = NULL ;

        //----------------------------------------------------------------------
        // U = triu(C)
        //----------------------------------------------------------------------

        GRB_TRY (GrB_select (U, NULL, NULL, GrB_TRIU, C, 0, NULL)) ;

    }
    else
    {
        GrB_Matrix AT = NULL, D = NULL ;

        GRB_TRY (GrB_Matrix_new (&AT, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_transpose (AT, NULL, GrB_ONEB_FP64, S, NULL)) ;

        //----------------------------------------------------------------------
        // C = A \/ A' to create an undirected graph C
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_eWiseAdd (C, NULL, NULL, GrB_ONEB_FP64, S, AT, NULL)) ;

        //----------------------------------------------------------------------
        // D = A + A' to create an undirected multigraph D
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Matrix_new (&D, GrB_FP64, n, n)) ;
        GRB_TRY (GrB_eWiseAdd (D, NULL, NULL, GrB_PLUS_FP64, S, AT, NULL)) ;

        GRB_TRY (GrB_free (&AT)) ;
        GRB_TRY (GrB_free (&S)) ;

        //----------------------------------------------------------------------
        // U = triu(D)
        //----------------------------------------------------------------------

        // note that L=U' since D is symmetric
        GRB_TRY (GrB_select (U, NULL, NULL, GrB_TRIU, D, 0, NULL)) ;
        GRB_TRY (GrB_free (&D)) ;
    }

    //--------------------------------------------------------------------------
    // Find wedges of each node
    //--------------------------------------------------------------------------

    // W(i) = sum (C (i,:)) = # of entries in C(i,:) because all entries
    // present in C are equal to 1.
    GRB_TRY (GrB_Vector_new (&W, GrB_FP64, n)) ;
    // x = zeros (n,1)
    GRB_TRY (GrB_Vector_new (&x, GrB_INT64, n)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    // W = C*x using the plus_one semiring
    GRB_TRY (GrB_mxv (W, NULL, NULL, LAGraph_plus_one_fp64, C, x, NULL)) ;
    GrB_free (&x) ;

    // Compute vector W defining the number of wedges per vertex
    GRB_TRY (GrB_apply(W, NULL, NULL, LAGraph_COMB_FP64, W, NULL));

    //--------------------------------------------------------------------------
    // Calculate triangles
    //--------------------------------------------------------------------------

    // CL<C> = C*L = C*U' using a masked dot product
    GRB_TRY (GrB_Matrix_new (&CL, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_mxm (CL, C, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C, U,
                         GrB_DESC_ST1));
    GRB_TRY (GrB_free (&U)) ;

    //--------------------------------------------------------------------------
    // Calculate LCC
    //--------------------------------------------------------------------------

    // LCC(i) = sum (CL (i,:)) = # of triangles at each node
    GRB_TRY (GrB_Vector_new (&LCC, GrB_FP64, n)) ;
    GRB_TRY (GrB_reduce (LCC, NULL, NULL, GrB_PLUS_FP64, CL, NULL)) ;
    GRB_TRY (GrB_free (&CL)) ;

    // LCC = LCC ./ W
    GRB_TRY (GrB_eWiseMult (LCC, NULL, NULL, GrB_DIV_FP64, LCC, W, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*LCC_handle) = LCC ; LCC = NULL ;

    LG_FREE_ALL ;
    return (GrB_SUCCESS) ;
#endif
}
