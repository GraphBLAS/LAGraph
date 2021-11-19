//------------------------------------------------------------------------------
// LAGraph_msf.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

/**
 * Code is based on Boruvka's minimum spanning forest algorithm
 * Contributed by Yongzhe Zhang (zyz915@gmail.com)
 */

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

#define LAGraph_FREE_ALL                             \
{                                                    \
    GrB_free (&S);                                   \
    GrB_free (&T);                                   \
    free(I); free(V);                                \
    free(SI); free(SJ); free(SX);                    \
    free(parent); free(partner); free(weight);       \
    GrB_free (&S);                      \
    GrB_free (&f);                      \
    GrB_free (&i);                      \
    GrB_free (&t);                      \
    GrB_free (&edge);                   \
    GrB_free (&cedge);                  \
    GrB_free (&mask);                   \
    GrB_free (&index);                  \
    GrB_free (&comb);                   \
    GrB_free (&combMin);                \
    GrB_free (&fst);                    \
    GrB_free (&snd);                    \
    GrB_free (&s1);                     \
    GrB_free (&s2);                     \
}

#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
// encode each edge into a single uint64_t
static void combine (void *z, const void *x, const void *y)
{
    *(uint64_t*)z = ((*(uint64_t*)x) << 32) + (*(uint64_t*)y);
}

static void get_fst (void *y, const void *x)
{
    *(uint64_t*)y = (*(uint64_t*)x) >> 32;
}

static void get_snd (void *y, const void *x)
{
    *(uint64_t*)y = (*(uint64_t*)x) & INT_MAX;
}

//****************************************************************************
// w[index[i]] = min(w[index[i]], s[i]) for i in [0..n-1]
static GrB_Info Reduce_assign (GrB_Vector w,
        GrB_Vector s, GrB_Index *index, GrB_Index n)
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
// global C arrays (for implementing various GxB_SelectOp)
static GrB_Index *weight = NULL, *parent = NULL, *partner = NULL;

// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component
bool f1 (
    GrB_Index i, GrB_Index j,
    const void *x, const void *thunk)
{
    uint64_t *aij = (uint64_t*) x;
    return (weight[i] == *aij) && (parent[j] == partner[i]);
}

// edge removal:
// A(i, j) is removed when parent[i] == parent[j]
bool f2 (
    GrB_Index i, GrB_Index j,
    const void *x, const void *thunk)
{
    uint64_t *aij = (uint64_t*) x;
    return (parent[i] != parent[j]);
}

//****************************************************************************
//****************************************************************************
GrB_Info LAGraph_msf
(
    GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,       // input matrix
    bool sanitize       // if true, ensure A is symmetric
)
{
#if !LG_SUITESPARSE
    // currently requires GxB_select; FIXME: make pure GrB with GrB_select
    return GrB_PANIC;
#else
    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, i = NULL, t = NULL,
        edge = NULL, cedge = NULL, mask = NULL, index = NULL;
    GrB_Index *I = NULL, *V = NULL, *SI = NULL, *SJ = NULL, *SX = NULL;

    GrB_BinaryOp comb = NULL;
    GrB_Semiring combMin = NULL;
    GrB_UnaryOp fst = NULL, snd = NULL;

    GxB_SelectOp s1 = NULL, s2 = NULL;

    LAGRAPH_OK (GrB_Matrix_nrows (&n, A));

    if (sanitize)
    {
        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_UINT64, n, n));
        LAGRAPH_OK (GrB_eWiseAdd (S, 0, 0, GrB_MIN_UINT64, A, A, GrB_DESC_T1));
    }
    else
    {
        // Use the input as-is, and assume it is GrB_UINT64 and symmetric
        LAGRAPH_OK (GrB_Matrix_dup (&S, A));
    }

    LAGRAPH_OK (GrB_Matrix_new (&T, GrB_UINT64, n, n));

    LAGRAPH_OK (GrB_Vector_new (&t, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&f, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&i, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&edge, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&cedge, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_BOOL, n));
    LAGRAPH_OK (GrB_Vector_new (&index, GrB_UINT64, n));

    // temporary arrays
    I = malloc (sizeof(GrB_Index) * n);
    V = malloc (sizeof(GrB_Index) * n);
    SI = malloc (sizeof(GrB_Index) * n * 2);
    SJ = malloc (sizeof(GrB_Index) * n * 2);
    SX = malloc (sizeof(GrB_Index) * n * 2);

    // global arrays
    parent = malloc (sizeof(GrB_Index) * n);
    weight = malloc (sizeof(GrB_Index) * n);
    partner = malloc (sizeof(GrB_Index) * n);

    // prepare vectors
    for (GrB_Index i = 0; i < n; i++)
        I[i] = parent[i] = i;
    LAGRAPH_OK (GrB_Vector_build (f, I, parent, n, GrB_PLUS_UINT64));
    LAGRAPH_OK (GrB_assign (i, 0, 0, f, GrB_ALL, 0, 0));

    // semiring & monoid
    GrB_Index inf = ((uint64_t) INT_MAX << 32) ^ INT_MAX;
    LAGRAPH_OK (GrB_BinaryOp_new (&comb, combine, GrB_UINT64, GrB_UINT64, GrB_UINT64));
    LAGRAPH_OK (GrB_Semiring_new (&combMin, GrB_MIN_MONOID_UINT64, comb));
    LAGRAPH_OK (GrB_UnaryOp_new (&fst, get_fst, GrB_UINT64, GrB_UINT64));
    LAGRAPH_OK (GrB_UnaryOp_new (&snd, get_snd, GrB_UINT64, GrB_UINT64));

    // GrB_SelectOp
    GxB_SelectOp_new (&s1, f1, GrB_UINT64, GrB_NULL);
    GxB_SelectOp_new (&s2, f2, GrB_UINT64, GrB_NULL);

    // the main computation
    GrB_Index nvals, diff, ntuples = 0, num;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, S));
    for (int iters = 1; nvals > 0; iters++)
    {
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        LAGRAPH_OK (GrB_assign (edge, 0, 0, inf, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_mxv (edge, 0, GrB_MIN_UINT64, combMin, S, f, 0));
        // cedge[u] = children's minimum edge  | if u is a root
        //          = (INT_MAX, u)             | otherwise
        LAGRAPH_OK (GrB_assign (t, 0, 0, (uint64_t) INT_MAX, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_eWiseMult (cedge, 0, 0, comb, t, i, 0));
        LAGRAPH_OK (Reduce_assign (cedge, edge, parent, n));
        // if (f[u] == u) f[u] := snd(cedge[u])  -- the index part of the edge
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, f, i, 0));
        LAGRAPH_OK (GrB_apply (f, mask, GrB_SECOND_UINT64, snd, cedge, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
        LAGRAPH_OK (GrB_extract (t, 0, 0, f, V, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, i, t, 0));
        LAGRAPH_OK (GrB_assign (f, mask, GrB_MIN_UINT64, i, GrB_ALL, 0, 0));

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, i, f, 0));
        LAGRAPH_OK (GrB_assign (cedge, mask, 0, inf, GrB_ALL, 0, 0));

        // 2. every vertex tries to know whether one of its edges is selected
        LAGRAPH_OK (GrB_extract (t, 0, 0, cedge, parent, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, edge, t, 0));

        // 3. each root picks a vertex from its children to generate the solution
        LAGRAPH_OK (GrB_assign (index, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (index, mask, 0, i, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (t, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (Reduce_assign (t, index, parent, n));
        LAGRAPH_OK (GrB_extract (index, 0, 0, t, parent, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, i, index, 0));

        // 4. generate the select function (set the global pointers)
        LAGRAPH_OK (GrB_assign (t, 0, 0, inf, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_apply (t, mask, 0, fst, edge, 0));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, weight, &n, t));
        LAGRAPH_OK (GrB_assign (t, 0, 0, inf, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_apply (t, mask, 0, snd, edge, 0));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, partner, &n, t));
        LAGRAPH_OK (GxB_select (T, 0, 0, s1, S, GrB_NULL, 0));
        LAGRAPH_OK (GrB_Vector_clear (t));

        // 5. the generated matrix may still have redundant edges
        //    remove the duplicates by GrB_mxv() and store them as tuples
        LAGRAPH_OK (GrB_Vector_clear (edge));
        LAGRAPH_OK (GrB_mxv (edge, mask, GrB_MIN_UINT64, combMin, T, i, 0));
        LAGRAPH_OK (GrB_Vector_nvals (&num, edge));
        LAGRAPH_OK (GrB_apply (t, 0, 0, snd, edge, 0));
        LAGRAPH_OK (GrB_Vector_extractTuples (SI + ntuples, SJ + ntuples, &num, t));
        LAGRAPH_OK (GrB_apply (t, 0, 0, fst, edge, 0));
        LAGRAPH_OK (GrB_Vector_extractTuples (SI + ntuples, SX + ntuples, &num, t));
        LAGRAPH_OK (GrB_Vector_clear (t));
        ntuples += num;

        // path halving until every vertex points on a root
        do {
            LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
            LAGRAPH_OK (GrB_extract (t, 0, 0, f, V, n, 0));
            LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, f, t, 0));
            LAGRAPH_OK (GrB_assign (f, 0, 0, t, GrB_ALL, 0, 0));
            LAGRAPH_OK (GrB_reduce (&diff, 0, GrB_PLUS_MONOID_UINT64, mask, 0));
        } while (diff != 0);

        // remove the edges in the same connected component
        LAGRAPH_OK (GrB_Vector_extractTuples (I, parent, &n, f));
        LAGRAPH_OK (GxB_select (S, 0, 0, s2, S, GrB_NULL, 0));
        GrB_Matrix_nvals (&nvals, S);
        if (nvals == 0) break;
    }
    LAGRAPH_OK (GrB_Matrix_clear (T));
    LAGRAPH_OK (GrB_Matrix_build (T, SI, SJ, SX, ntuples, GrB_SECOND_UINT64));
    *result = T;

    LAGraph_FREE_ALL;
    return GrB_SUCCESS;
#endif
}
