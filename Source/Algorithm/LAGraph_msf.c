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

/**
 * Code is based on Boruvka's minimum spanning forest algorithm
 * Contributed by Yongzhe Zhang (zyz915@gmail.com)
 */

#define LAGRAPH_FREE_ALL

#include "LAGraph.h"

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

// global C arrays (for implementing various GxB_SelectOp)
static GrB_Index *weight, *parent, *partner;
// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component
static bool f1 (GrB_Index i, GrB_Index j, GrB_Index nrows, GrB_Index ncols,
        const void *x, const void *thunk)
{
    uint64_t *aij = (uint64_t*) x;
    return (weight[i] == *aij) && (parent[j] == partner[i]);
}
// edge removal:
// A(i, j) is removed when parent[i] == parent[j]
static bool f2 (GrB_Index i, GrB_Index j, GrB_Index nrows, GrB_Index ncols,
        const void *x, const void *thunk)
{
    uint64_t *aij = (uint64_t*) x;
    return (parent[i] != parent[j]);
}

GrB_Info LAGraph_msf
(
    GrB_Matrix *result,     // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info;
    GrB_Index n;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A));

    GrB_Matrix S;
    if (sanitize) {
        GrB_Descriptor desc;
        LAGRAPH_OK (GrB_Descriptor_new (&desc));
        LAGRAPH_OK (GrB_Descriptor_set (desc, GrB_INP1, GrB_TRAN));
        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_UINT64, n, n));
        LAGRAPH_OK (GrB_eWiseAdd (S, 0, 0, GrB_MIN_UINT64, A, A, desc));
        LAGRAPH_OK (GrB_free (&desc));
    } else {
        // Use the input as-is, and assume it is GrB_UINT64 and symmetric
        LAGRAPH_OK (GrB_Matrix_dup (&S, A));
    }
    // matrix
    GrB_Matrix T;
    LAGRAPH_OK (GrB_Matrix_new (&T, GrB_UINT64, n, n));
    // vectors
    GrB_Vector f, i, t, edge, cedge, mask, index;
    LAGRAPH_OK (GrB_Vector_new (&t, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&f, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&i, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&edge, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&cedge, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&index, GrB_UINT64, n));
    // temporary arrays
    GrB_Index *I = malloc (sizeof(GrB_Index) * n);
    GrB_Index *V = malloc (sizeof(GrB_Index) * n);
    GrB_Index *SI = malloc (sizeof(GrB_Index) * n * 2);
    GrB_Index *SJ = malloc (sizeof(GrB_Index) * n * 2);
    GrB_Index *SX = malloc (sizeof(GrB_Index) * n * 2);
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
    GrB_Monoid Min, Add;
    GrB_BinaryOp comb;
    GrB_Semiring combMin;
    GrB_UnaryOp fst, snd;
    GrB_Index inf = ((uint64_t) INT_MAX << 32) ^ INT_MAX;
    LAGRAPH_OK (GrB_BinaryOp_new (&comb, combine, GrB_UINT64, GrB_UINT64, GrB_UINT64));
    LAGRAPH_OK (GrB_Monoid_new (&Min, GrB_MIN_UINT64, inf));
    LAGRAPH_OK (GrB_Monoid_new (&Add, GrB_PLUS_UINT64, (uint64_t)0));
    LAGRAPH_OK (GrB_Semiring_new (&combMin, Min, comb));
    LAGRAPH_OK (GrB_UnaryOp_new (&fst, get_fst, GrB_UINT64, GrB_UINT64));
    LAGRAPH_OK (GrB_UnaryOp_new (&snd, get_snd, GrB_UINT64, GrB_UINT64));
    // GrB_SelectOp
    GxB_SelectOp s1, s2;
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
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, f, i, 0));
        LAGRAPH_OK (GrB_apply (f, mask, GrB_SECOND_UINT64, snd, cedge, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
        LAGRAPH_OK (GrB_extract (t, 0, 0, f, V, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, i, t, 0)); 
        LAGRAPH_OK (GrB_assign (f, mask, GrB_MIN_UINT64, i, GrB_ALL, 0, 0));
        // five steps to generate the solution 
        // 1. new roots (f[i] == i) revise their entries in cedge
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, i, f, 0));
        LAGRAPH_OK (GrB_assign (cedge, mask, 0, inf, GrB_ALL, 0, 0));
        // 2. every vertex tries to know whether one of its edges is selected
        LAGRAPH_OK (GrB_extract (t, 0, 0, cedge, parent, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask ,0, 0, GxB_ISEQ_UINT64, edge, t, 0));
        // 3. each root picks a vertex from its children to generate the solution
        LAGRAPH_OK (GrB_assign (index, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (index, mask, 0, i, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (t, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (Reduce_assign (t, index, parent, n));
        LAGRAPH_OK (GrB_extract (index, 0, 0, t, parent, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask ,0, 0, GxB_ISEQ_UINT64, i, index, 0));
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
            LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISNE_UINT64, f, t, 0));
            LAGRAPH_OK (GrB_assign (f, 0, 0, t, GrB_ALL, 0, 0));
            LAGRAPH_OK (GrB_reduce (&diff, 0, Add, mask, 0));
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
    // free
    free(I); free(V);
    free(SI); free(SJ); free(SX);
    free(parent); free(partner); free(weight);

    LAGRAPH_OK (GrB_free (&S));
    LAGRAPH_OK (GrB_free (&f));
    LAGRAPH_OK (GrB_free (&i));
    LAGRAPH_OK (GrB_free (&t));
    LAGRAPH_OK (GrB_free (&edge));
    LAGRAPH_OK (GrB_free (&cedge));
    LAGRAPH_OK (GrB_free (&mask));
    LAGRAPH_OK (GrB_free (&index));
    LAGRAPH_OK (GrB_free (&Min));
    LAGRAPH_OK (GrB_free (&Add));
    LAGRAPH_OK (GrB_free (&comb));
    LAGRAPH_OK (GrB_free (&combMin));
    LAGRAPH_OK (GrB_free (&fst));
    LAGRAPH_OK (GrB_free (&snd));
    LAGRAPH_OK (GrB_free (&s1));
    LAGRAPH_OK (GrB_free (&s2));
    return GrB_SUCCESS;
}
