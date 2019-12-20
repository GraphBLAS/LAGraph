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
 * Code is based on Borilvka's minimum spanning forest algorithm
 * Author: Yongzhe Zhang
 **/

#define LAGRAPH_FREE_ALL

#include "LAGraph.h"

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

static GrB_Index *I, *V;
static bool func(const GrB_Index i, const GrB_Index j, const GrB_Index nrows,
        const GrB_Index ncols, const void *x, const void *thunk)
{
    return V[i] != V[j];
}

GrB_Info LAGraph_cc_boruvka
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info;
    GrB_Index n;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    GrB_Matrix S;
    if (sanitize)
    {
        GrB_Descriptor desc;
        LAGRAPH_OK (GrB_Descriptor_new (&desc));
        LAGRAPH_OK (GrB_Descriptor_set (desc, GrB_INP1, GrB_TRAN));
        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_BOOL, n, n));
        LAGRAPH_OK (GrB_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, desc));
        LAGRAPH_FREE (desc);
    }
    else
    {
        // Use the input as-is, and assume it is binary and symmetric
        LAGRAPH_OK (GrB_Matrix_dup (&S, A));
    }
    GrB_Vector f, gp, mnp, ccmn, i, inf, mask;
    LAGRAPH_OK (GrB_Vector_new (&f, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&gp, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mnp, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&ccmn, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&i, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&inf, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_UINT64, n));
    // prepare
    I = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    V = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    for (GrB_Index i = 0; i < n; i++)
        I[i] = V[i] = i;
    LAGRAPH_OK (GrB_Vector_build (f, I, V, n, GrB_PLUS_UINT64));
    LAGRAPH_OK (GrB_assign (i, 0, 0, f, GrB_ALL, 0, 0));
    LAGRAPH_OK (GrB_assign (inf, 0, 0, n, GrB_ALL, 0, 0));
    // semiring & monoid
    GrB_Monoid Min, Add;
    GrB_Semiring sel2ndMin;
    LAGRAPH_OK (GrB_Monoid_new (&Min, GrB_MIN_UINT64, n));
    LAGRAPH_OK (GrB_Semiring_new (&sel2ndMin, Min, GrB_SECOND_UINT64));
    LAGRAPH_OK (GrB_Monoid_new (&Add, GrB_PLUS_UINT64, (GrB_Index) 0));
    GrB_Index nvals, diff;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, A));
    // GxB_SelectOp
    GxB_SelectOp sel_op;
    LAGRAPH_OK (GxB_SelectOp_new (&sel_op, func, 0, 0));
    for (int iters = 1; nvals > 0; iters++)
    {
        // every vertex points t oa root vertex at the begining
        // mnp[u] = u's minimum neighbor's parent
        LAGRAPH_OK (GrB_assign (mnp, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_mxv (mnp, 0, GrB_MIN_UINT64, sel2ndMin, S, f, 0));
        // ccmn[u] = connect component's minimum neighbor | if u is a root
        //         = inf                                  | otherwise
        LAGRAPH_OK (GrB_assign (ccmn, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (Reduce_assign (ccmn, mnp, V, n));
        // f[u] = ccmn[u] if ccmn[u] != inf
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISNE_UINT64, ccmn, inf, 0));
        LAGRAPH_OK (GrB_assign (f, mask, 0, ccmn, GrB_ALL, 0, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
        LAGRAPH_OK (GrB_extract (gp, 0, 0, f, V, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISEQ_UINT64, i, gp, 0)); 
        LAGRAPH_OK (GrB_assign (f, mask, GrB_MIN_UINT64, i, GrB_ALL, 0, 0));
        // shortcutting f[i] = f[f[i]]
        do {
            LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, f));
            LAGRAPH_OK (GrB_extract (gp, 0, 0, f, V, n, 0));
            LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GxB_ISNE_UINT64, f, gp, 0));
            LAGRAPH_OK (GrB_assign (f, 0, 0, gp, GrB_ALL, 0, 0));
            LAGRAPH_OK (GrB_reduce (&diff, 0, Add, mask, 0));
        } while (diff != 0);
        // remove edges connecting the same connected component
        LAGRAPH_OK (GxB_select (S, 0, 0, sel_op, S, 0, 0));
        LAGRAPH_OK (GrB_Matrix_nvals (&nvals, S));
    }
    *result = f;

    free(I);
    free(V);
    GrB_free (&gp);
    GrB_free (&mnp);
    GrB_free (&ccmn);
    GrB_free (&i);
    GrB_free (&inf);
    GrB_free (&mask);
    GrB_free (&Add);
    GrB_free (&Min);
    GrB_free (&sel2ndMin);
    return GrB_SUCCESS;
}
