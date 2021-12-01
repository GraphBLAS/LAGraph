//------------------------------------------------------------------------------
// LAGraph_scc.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

/**
 * Code is based on the Min-Label algorithm described in the following paper:
 * D. Yan, J. Cheng, K. Xin, Y. Lu, W. Ng, Y. Bu, "Pregel Algorithms for Graph
 * Connectivity Problems with Performance Guarantees"
 * Proc. VLDB Endow. 7, 14 (October 2014), 1821–1832.
 * DOI: https://doi.org/10.14778/2733085.2733089
 *
 * Contributed by Yongzhe Zhang (zyz915@gmail.com)
 **/

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

#define LAGraph_FREE_ALL ;

#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
// global C arrays used in SelectOp
GrB_Index *I = NULL, *V = NULL, *F = NULL, *B = NULL, *M = NULL;

// edge_removal:
//  - remove the edges connected to newly identified SCCs (vertices u with M[u]==1)
//  - remove the edges (u, v) where u and v can never be in the same SCC.
//
// Here's a brief explanation of the second case. After the forward and backward
// propagation, each vertex u has two labels
//  - F[u]: the smallest vertex that can reach u
//  - B[u]: the smallest vertex that is reachable from u
// If two vertices u and v are in the same SCC, then F[u]==F[v] and B[u]==B[v] must
// hold. The converse is not true unless F[u]==B[u]. However, we can safely remove
// an edge (u, v) if either F[u]!=F[v] or B[u]!=B[v] holds, which can accelerate
// the SCC computation in the future rounds.
bool edge_removal (GrB_Index i, GrB_Index j, const void *x, const void *thunk) ;
bool edge_removal (GrB_Index i, GrB_Index j, const void *x, const void *thunk)
{
    return !M[i] && !M[j] && F[i] == F[j] && B[i] == B[j];
}

//****************************************************************************
// trim_one: remove the edges connected to trivial SCCs
//  - A vertex is a trivial SCC if it has no incoming or outgoing edges.
//  - M[i] = i   | if vertex i is a trivial SCC
//    M[i] = n   | otherwise
bool trim_one (GrB_Index i, GrB_Index j, const void *x, const void *thunk) ;
bool trim_one (GrB_Index i, GrB_Index j, const void *x, const void *thunk)
{
    return M[i] == M[j];
}

//****************************************************************************
// label propagation
//  - label  : (input/output) labels
//  - mask   : (input) mask
//  - A      : (input) original matrix
//  - AT     : (input) transposed matrix
//  - n      : (input) number of vertices
//  - is_csr : (input) format
static GrB_Info propagate (GrB_Vector label, GrB_Vector mask,
        GrB_Matrix A, GrB_Matrix AT, GrB_Index n, bool is_csr)
{
    GrB_Info info;
    // semirings

    GrB_Vector s, t;
    LAGRAPH_OK (GrB_Vector_new (&s, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&t, GrB_UINT64, n));
    LAGRAPH_OK (GrB_assign (s, mask, 0, label, GrB_ALL, 0, 0));
    LAGRAPH_OK (GrB_assign (t, 0, 0, label, GrB_ALL, 0, 0));

    GrB_Index active;
    while (true)
    {
        if (is_csr) {
            LAGRAPH_OK (GrB_vxm (t, 0, GrB_MIN_UINT64,
                                 GrB_MIN_FIRST_SEMIRING_UINT64, s, A, 0));
        } else {
            LAGRAPH_OK (GrB_mxv (t, 0, GrB_MIN_UINT64,
                                 GrB_MIN_SECOND_SEMIRING_UINT64, AT, s, 0));
        }
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISNE_UINT64, t, label, 0));
        LAGRAPH_OK (GrB_assign (label, mask, 0, t, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_reduce (&active, 0, GrB_PLUS_MONOID_UINT64, mask, 0));
        if (active == 0) break;
        LAGRAPH_OK (GrB_Vector_clear (s));
        LAGRAPH_OK (GrB_assign (s, mask, 0, label, GrB_ALL, 0, 0));
    }

    LAGRAPH_OK (GrB_free (&s));
    LAGRAPH_OK (GrB_free (&t));
    return GrB_SUCCESS;
}

//****************************************************************************
//****************************************************************************
GrB_Info LAGraph_scc
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A            // input matrix
)
{
#if !(LG_SUITESPARSE)
    // FIXME: can be GrB*
    return GrB_PANIC;
#else
    GrB_Info info;
    GrB_Vector scc;
    GrB_Vector ind;
    GrB_Vector inf;
    GrB_Vector f, b, mask;
    GxB_SelectOp sel1, sel2;
    GrB_Monoid Add;

    GrB_Index n, nvals;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A));

    // store the graph in both directions (forward / backward)
    GrB_Matrix FW, BW;
    LAGRAPH_OK (GrB_Matrix_new (&FW, GrB_BOOL, n, n));
    LAGRAPH_OK (GrB_Matrix_new (&BW, GrB_BOOL, n, n));
    LAGRAPH_OK (GrB_transpose (FW, 0, 0, A, GrB_DESC_T0)); // FW = A
    LAGRAPH_OK (GrB_transpose (BW, 0, 0, A, 0));     // BW = A'

    // check format
    GxB_Format_Value A_format, AT_format;
    LAGRAPH_OK (GxB_get (FW, GxB_FORMAT, &A_format));
    LAGRAPH_OK (GxB_get (BW, GxB_FORMAT, &AT_format));

    bool is_csr = (A_format == GxB_BY_ROW && AT_format == GxB_BY_ROW);
    bool is_csc = (A_format == GxB_BY_COL && AT_format == GxB_BY_COL);
    if (!is_csr && !is_csc) {
        LAGRAPH_ERROR ("A and AT must in the same format:\n"
            "both GxB_BY_ROW, or both GxB_BY_COL",
            GrB_INVALID_VALUE) ;
    }

    I = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    V = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    F = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    B = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    M = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    for (GrB_Index i = 0; i < n; i++)
        I[i] = V[i] = i;

    // scc: the SCC identifier for each vertex
    // scc[u] == n: not assigned yet
    LAGRAPH_OK (GrB_Vector_new (&scc, GrB_UINT64, n));
    // vector of indices: ind[i] == i
    LAGRAPH_OK (GrB_Vector_new (&ind, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_build (ind, I, V, n, GrB_PLUS_UINT64));
    // vector of infinite value: inf[i] == n
    LAGRAPH_OK (GrB_Vector_new (&inf, GrB_UINT64, n));
    LAGRAPH_OK (GrB_assign (inf, 0, 0, n, GrB_ALL, 0, 0));
    // other vectors
    LAGRAPH_OK (GrB_Vector_new (&f, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&b, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_UINT64, n));
    // GxB_SelectOp
    LAGRAPH_OK (GxB_SelectOp_new (&sel1, trim_one, GrB_BOOL, GrB_NULL));
    LAGRAPH_OK (GxB_SelectOp_new (&sel2, edge_removal, GrB_BOOL, GrB_NULL));

    // remove trivial SCCs
    LAGRAPH_OK (GrB_reduce (f, 0, GrB_PLUS_UINT64, GrB_PLUS_UINT64, FW, 0));
    LAGRAPH_OK (GrB_reduce (b, 0, GrB_PLUS_UINT64, GrB_PLUS_UINT64, BW, 0));
    LAGRAPH_OK (GrB_eWiseMult (mask, 0, GxB_LAND_UINT64, GxB_LAND_UINT64, f, b, 0));
    LAGRAPH_OK (GrB_Vector_nvals (&nvals, mask));

    LAGRAPH_OK (GrB_assign (scc, 0, 0, ind, GrB_ALL, 0, 0));
    LAGRAPH_OK (GrB_assign (scc, mask, 0, n, GrB_ALL, 0, 0));
    LAGRAPH_OK (GrB_Vector_clear (mask));

    if (nvals < n)
    {
        LAGRAPH_OK (GrB_Vector_extractTuples (I, M, &n, scc));
        LAGRAPH_OK (GxB_select (FW, 0, 0, sel1, FW, GrB_NULL, 0));
        LAGRAPH_OK (GxB_select (BW, 0, 0, sel1, BW, GrB_NULL, 0));
    }

    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, FW));
    while (nvals > 0)
    {
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, scc, inf, 0));
        LAGRAPH_OK (GrB_assign (f, 0, 0, ind, GrB_ALL, 0, 0));
        LAGRAPH_OK (propagate (f, mask, FW, BW, n, is_csr));

        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, f, ind, 0));
        LAGRAPH_OK (GrB_assign (b, 0, 0, inf, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (b, mask, 0, ind, GrB_ALL, 0, 0));
        LAGRAPH_OK (propagate (b, mask, BW, FW, n, is_csr));

        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, f, b, 0));
        LAGRAPH_OK (GrB_assign (scc, mask, GrB_MIN_UINT64, f, GrB_ALL, 0, 0));

        LAGRAPH_OK (GrB_Vector_extractTuples (I, F, &n, f));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, B, &n, b));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, M, &n, mask));

        LAGRAPH_OK (GxB_select (FW, 0, 0, sel2, FW, GrB_NULL, 0));
        LAGRAPH_OK (GxB_select (BW, 0, 0, sel2, BW, GrB_NULL, 0));

        LAGRAPH_OK (GrB_Matrix_nvals (&nvals, FW));
    }
    LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, scc, inf, 0));
    LAGRAPH_OK (GrB_assign (scc, mask, 0, ind, GrB_ALL, 0, 0));

    LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, scc, ind, 0));
    LAGRAPH_OK (GrB_reduce (&nvals, 0, GrB_PLUS_MONOID_UINT64, mask, 0));

    *result = scc;
    scc = NULL;

    free (I);
    free (V);
    free (F);
    free (B);
    free (M);
    GrB_free (&ind);
    GrB_free (&inf);
    GrB_free (&f);
    GrB_free (&b);
    GrB_free (&mask);
    GrB_free (&FW);
    GrB_free (&BW);
    GrB_free (&sel1);
    GrB_free (&sel2);
    GrB_free (&scc);

    return GrB_SUCCESS;
#endif
}
