//------------------------------------------------------------------------------
// LAGraph_cc_lacc.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

/**
 * Code is based on the algorithm described in the following paper
 * Azad, Buluç LACC: a linear-algebraic algorithm for finding connected components
 * in distributed memory (IPDPS 2019)
 **/

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL                \
    free(I);                            \
    free(V);                            \
    GrB_free (&stars);     \
    GrB_free (&mask);      \
    GrB_free (&parents);   \
    GrB_free (&gp);        \
    GrB_free (&mnp);       \
    GrB_free (&hookMNP);   \
    GrB_free (&hookP);     \
    GrB_free (&pNonstars); \
    GrB_free (&tmp);       \
    GrB_free (&nsgp);

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

//****************************************************************************
// mask = NULL, accumulator = GrB_MIN_UINT64, descriptor = NULL
static GrB_Info Reduce_assign (GrB_Vector w,
                               GrB_Vector src,
                               GrB_Index *index,
                               GrB_Index nLocs)
{
    GrB_Index nw, ns;
    GrB_Vector_nvals(&nw, w);
    GrB_Vector_nvals(&ns, src);
    GrB_Index *mem = (GrB_Index*) malloc(sizeof(GrB_Index) * nw * 3);
    GrB_Index *ind = mem, *sval = mem + nw, *wval = sval + nw;
    GrB_Vector_extractTuples(ind, wval, &nw, w);
    GrB_Vector_extractTuples(ind, sval, &ns, src);
    for (GrB_Index i = 0; i < nLocs; i++)
        if (sval[i] < wval[index[i]])
            wval[index[i]] = sval[i];
    GrB_Vector_clear(w);
    GrB_Vector_build(w, ind, wval, nw, GrB_PLUS_UINT64);
    free(mem);
    return GrB_SUCCESS;
}

//****************************************************************************
GrB_Info LAGraph_cc_lacc
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    if (result == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    GrB_Info info;

    GrB_Vector stars = NULL, mask = NULL;
    GrB_Vector parents = NULL, gp = NULL, mnp = NULL;
    GrB_Vector hookMNP = NULL, hookP = NULL;
    GrB_Vector tmp = NULL, pNonstars = NULL, nsgp = NULL; // temporary
    GrB_Index *I = NULL;
    GrB_Index *V = NULL;

    GrB_Index n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    //GrB_Index nnz ;
    //LAGRAPH_OK (GrB_Matrix_nvals (&nnz, A)) ;
    //printf ("number of nodes: %g\n", (double) n) ;
    //printf ("number of edges: %g\n", (double) nnz) ;

    GrB_Matrix S = NULL;
    if (sanitize)
    {
        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_BOOL, n, n)) ;
        LAGRAPH_OK (GrB_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, GrB_DESC_T1)) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary and symmetric
        S = A ;
    }

    // vectors
    LAGRAPH_OK (GrB_Vector_new (&stars, GrB_BOOL, n));
    LAGRAPH_OK (GrB_Vector_new (&mask, GrB_BOOL, n));
    LAGRAPH_OK (GrB_Vector_new (&parents, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&gp, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&hookMNP, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&hookP, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&pNonstars, GrB_UINT64, n));

    // temporary arrays
    I = malloc(sizeof(GrB_Index) * n);
    V = malloc(sizeof(GrB_Index) * n);

    // prepare the vectors
    for (GrB_Index i = 0 ; i < n ; i++)
        I[i] = V[i] = i;
    LAGRAPH_OK (GrB_Vector_build (parents, I, V, n, GrB_PLUS_UINT64));
    LAGRAPH_OK (GrB_Vector_dup (&mnp, parents));
    LAGRAPH_OK (GrB_assign (stars, 0, 0, true, GrB_ALL, 0, 0)) ;

    // main computation
    GrB_Index nHooks, nStars, nNonstars;
    while (true) {
        // ---------------------------------------------------------
        // CondHook(A, parents, stars);
        // ---------------------------------------------------------
        LAGRAPH_OK (GrB_mxv (mnp, 0, 0, GrB_MIN_SECOND_SEMIRING_UINT64,
                             S, parents, 0));
        LAGRAPH_OK (GrB_Vector_clear (mask));
        LAGRAPH_OK (GrB_eWiseMult(mask, stars, 0, GrB_LT_UINT64, mnp, parents, 0));
        LAGRAPH_OK (GrB_assign (hookMNP, mask, 0, mnp, GrB_ALL, n, 0));
        LAGRAPH_OK (GrB_eWiseMult (hookP, 0, 0, GrB_SECOND_UINT64, hookMNP, parents, 0));
        LAGRAPH_OK (GrB_Vector_clear (mnp));
        LAGRAPH_OK (GrB_Vector_nvals (&nHooks, hookP));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &nHooks, hookP));
        LAGRAPH_OK (GrB_Vector_new (&tmp, GrB_UINT64, nHooks));
        LAGRAPH_OK (GrB_extract (tmp, 0, 0, hookMNP, I, nHooks, 0));
        LAGRAPH_OK (Reduce_assign (parents, tmp, V, nHooks));
        LAGRAPH_OK (GrB_Vector_clear (tmp));
        // modify the stars vector
        LAGRAPH_OK (GrB_assign (stars, 0, 0, false, V, nHooks, 0));
        LAGRAPH_OK (GrB_extract (tmp, 0, 0, parents, V, nHooks, 0)); // extract modified parents
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &nHooks, tmp));
        LAGRAPH_OK (GrB_assign (stars, 0, 0, false, V, nHooks, 0));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, parents));
        LAGRAPH_OK (GrB_extract (mask, 0, 0, stars, V, n, 0));
        LAGRAPH_OK (GrB_assign (stars, 0, GrB_LAND, mask, GrB_ALL, 0, 0));
        // clean up
        LAGRAPH_OK (GrB_Vector_clear (hookMNP));
        LAGRAPH_OK (GrB_Vector_clear (hookP));
        LAGRAPH_OK (GrB_free (&tmp));
        // ---------------------------------------------------------
        // UnCondHook(A, parents, stars);
        // ---------------------------------------------------------
        LAGRAPH_OK (GrB_assign (pNonstars, 0, 0, parents, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (pNonstars, stars, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_mxv (hookMNP, stars, 0, GrB_MIN_SECOND_SEMIRING_UINT64,
                             S, pNonstars, 0));
        // select the valid elemenets (<n) of hookMNP
        LAGRAPH_OK (GrB_assign (pNonstars, 0, 0, n, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_LT_UINT64, hookMNP, pNonstars, 0));
        LAGRAPH_OK (GrB_eWiseMult (hookP, mask, 0, GrB_SECOND_UINT64, hookMNP, parents, 0));
        LAGRAPH_OK (GrB_Vector_nvals (&nHooks, hookP));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &nHooks, hookP));
        LAGRAPH_OK (GrB_Vector_new (&tmp, GrB_UINT64, nHooks));
        LAGRAPH_OK (GrB_extract (tmp, 0, 0, hookMNP, I, nHooks, 0));
        LAGRAPH_OK (GrB_assign (parents, 0, 0, n, V, nHooks, 0)); // !!
        LAGRAPH_OK (Reduce_assign (parents, tmp, V, nHooks));
        // modify the star vector
        LAGRAPH_OK (GrB_assign (stars, 0, 0, false, V, nHooks, 0));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, parents));
        LAGRAPH_OK (GrB_extract (mask, 0, 0, stars, V, n, 0));
        LAGRAPH_OK (GrB_assign (stars, 0, GrB_LAND, mask, GrB_ALL, 0, 0));
        // check termination
        LAGRAPH_OK (GrB_reduce (&nStars, 0, GrB_PLUS_MONOID_UINT64, stars, 0));
        if (nStars == n) break;
        // clean up
        LAGRAPH_OK (GrB_Vector_clear(hookMNP));
        LAGRAPH_OK (GrB_Vector_clear(hookP));
        LAGRAPH_OK (GrB_Vector_clear(pNonstars));
        LAGRAPH_OK (GrB_free (&tmp));
        // ---------------------------------------------------------
        // Shortcut(parents);
        // ---------------------------------------------------------
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, parents));
        LAGRAPH_OK (GrB_extract (gp, 0, 0, parents, V, n, 0));
        LAGRAPH_OK (GrB_assign (parents, 0, 0, gp, GrB_ALL, 0, 0));
        // ---------------------------------------------------------
        // StarCheck(parents, stars);
        // ---------------------------------------------------------
        // calculate grandparents
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, parents));
        LAGRAPH_OK (GrB_extract (gp, 0, 0, parents, V, n, 0));
        // identify vertices whose parent and grandparent are different
        LAGRAPH_OK (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, gp, parents, 0));
        LAGRAPH_OK (GrB_Vector_new (&nsgp, GrB_UINT64, n));
        LAGRAPH_OK (GrB_assign (nsgp, mask, 0, gp, GrB_ALL, 0, 0));
        // extract indices and values for assign
        LAGRAPH_OK (GrB_Vector_nvals (&nNonstars, nsgp));
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &nNonstars, nsgp));
        LAGRAPH_OK (GrB_free (&nsgp));
        LAGRAPH_OK (GrB_assign (stars, 0, 0, true, GrB_ALL, 0, 0));
        LAGRAPH_OK (GrB_assign (stars, 0, 0, false, I, nNonstars, 0));
        LAGRAPH_OK (GrB_assign (stars, 0, 0, false, V, nNonstars, 0));
        // extract indices and values for assign
        LAGRAPH_OK (GrB_Vector_extractTuples (I, V, &n, parents));
        LAGRAPH_OK (GrB_extract (mask, 0, 0, stars, V, n, 0));
        LAGRAPH_OK (GrB_assign (stars, 0, GrB_LAND, mask, GrB_ALL, 0, 0));
    }
    *result = parents;
    parents = NULL ;        // return parents (set to NULL so it isn't freed)

    LAGraph_FREE_ALL;
    return GrB_SUCCESS;
}
