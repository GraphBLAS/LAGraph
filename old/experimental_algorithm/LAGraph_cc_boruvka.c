//------------------------------------------------------------------------------
// LAGraph_cc_boruvka.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

/**
 * Code is based on Borilvka's minimum spanning forest algorithm
 * Contributed by Yongzhe Zhang (zyz915@gmail.com)
 **/

#include <LAGraph.h>
#include <LAGraphX.h>
#include "LG_internal.h"

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL      \
    free(I);                  \
    free(V);                  \
    GrB_free (&gp);        \
    GrB_free (&mnp);       \
    GrB_free (&ccmn);      \
    GrB_free (&i);         \
    GrB_free (&inf);       \
    GrB_free (&mask);      \
    GrB_free (&select_op);

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

//****************************************************************************
// w[index[i]] = min(w[index[i]], s[i]) for i in [0..n-1]
static GrB_Info Reduce_assign (GrB_Vector w,
                               GrB_Vector s,
                               GrB_Index *index,
                               GrB_Index n)
{
    GrB_Index *mem = (GrB_Index*) malloc(sizeof(GrB_Index) * n * 3);
    GrB_Index *ind = mem, *sval = mem + n, *wval = sval + n;
    GrB_Vector_extractTuples(ind, wval, &n, w);
    GrB_Vector_extractTuples(ind, sval, &n, s);
    for (GrB_Index i = 0; i < n; i++)
        if (sval[i] < wval[index[i]])
            wval[index[i]] = sval[i];
    GrB_Vector_clear(w);
    GrB_Vector_build(w, ind, wval, n, GrB_PLUS_UINT64);
    free(mem);
    return GrB_SUCCESS;
}

//****************************************************************************
static GrB_Index *I, *V;

bool select_func(
                 const GrB_Index i, const GrB_Index j,
                 const void *x, const void *thunk)
{
    return V[i] != V[j];
}

//****************************************************************************
GrB_Info LAGraph_cc_boruvka
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info;

    // TODO remove once GrB_select is available and GxB_select is ported
#if !LG_SUITESPARSE
    return GrB_PANIC;
#else
    GrB_Index n;
    GrB_Vector f = NULL, gp = NULL, mnp = NULL, ccmn = NULL, i = NULL, inf = NULL, mask = NULL;
    GxB_SelectOp select_op = NULL;

    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    GrB_Matrix S;
    if (sanitize)
    {
        GrB_Descriptor desc;
        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_BOOL, n, n));
        LAGRAPH_OK (GrB_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, GrB_DESC_T1));
    }
    else
    {
        // Use the input as-is, and assume it is binary and symmetric
        LAGRAPH_OK (GrB_Matrix_dup (&S, A));
    }

    LAGRAPH_OK (GrB_Vector_new (&f, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&gp, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mnp, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&ccmn, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&i, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&inf, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_BOOL, n));

    // prepare
    I = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    V = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    for (GrB_Index i = 0; i < n; i++)
        I[i] = V[i] = i;
    LAGRAPH_OK (GrB_Vector_build (f, I, V, n, GrB_PLUS_UINT64));
    LAGRAPH_OK (GrB_assign (i, 0, 0, f, GrB_ALL, 0, 0));
    LAGRAPH_OK (GrB_assign (inf, 0, 0, n, GrB_ALL, 0, 0));

    // GxB_SelectOp
    LAGRAPH_OK (GxB_SelectOp_new (&select_op, select_func, 0, 0));

    GrB_Index nvals, diff;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, S));
    for (int iters = 1; nvals > 0; iters++)
    {
        // every vertex points to a root vertex at the begining
        // mnp[u] = u's minimum neighbor's parent
        LAGRAPH_OK (GrB_assign (mnp, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_mxv (mnp, 0, GrB_MIN_UINT64,
                             GrB_MIN_SECOND_SEMIRING_UINT64, S, f, 0));

        // ccmn[u] = connect component's minimum neighbor | if u is a root
        //         = inf                                  | otherwise
        LAGRAPH_OK (GrB_assign (ccmn, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (Reduce_assign (ccmn, mnp, V, n));

        // f[u] = ccmn[u] if ccmn[u] != inf
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, ccmn, inf, 0));
        LAGRAPH_OK (GrB_assign (f, mask, 0, ccmn, GrB_ALL, 0, 0));

        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
        LAGRAPH_OK (GrB_extract (gp, 0, 0, f, V, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, i, gp, 0));
        LAGRAPH_OK (GrB_assign (f, mask, GrB_MIN_UINT64, i, GrB_ALL, 0, 0));

        // shortcutting f[i] = f[f[i]]
        do {
            LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
            LAGRAPH_OK (GrB_extract (gp, 0, 0, f, V, n, 0));
            LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISNE_UINT64, f, gp, 0));
            LAGRAPH_OK (GrB_assign (f, 0, 0, gp, GrB_ALL, 0, 0));
            LAGRAPH_OK (GrB_reduce (&diff, 0, GrB_PLUS_MONOID_UINT64, mask, 0));
        } while (diff != 0);

        // remove the edges inside each connected component
        LAGRAPH_OK (GxB_select (S, 0, 0, select_op, S, 0, 0));
        LAGRAPH_OK (GrB_Matrix_nvals (&nvals, S));
    }
    *result = f;
    LAGraph_FREE_ALL;
    return GrB_SUCCESS;
#endif
}
