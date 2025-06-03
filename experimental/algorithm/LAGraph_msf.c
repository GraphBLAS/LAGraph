//------------------------------------------------------------------------------
// LAGraph_msf.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Yongzhe Zhang (zyz915@gmail.com)

//------------------------------------------------------------------------------

/**
 * Code is based on Boruvka's minimum spanning forest algorithm
 */

// TODO: is this ready for src?
// TODO: Reduce_assign is slow.  See src/algorithm/LG_CC_FastSV6/7.

#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
// encode each edge into a single uint64_t
static void combine (uint64_t *z, const uint64_t *x, const uint64_t *y)
{
    *z = ((*x) << 32) + (*y);
}

static void get_fst (uint64_t *y, const uint64_t *x)
{
    *y = (*x) >> 32;
}

static void get_snd (uint64_t *y, const uint64_t *x)
{
    *y = (*x) & UINT32_MAX;
}

//****************************************************************************
#undef  LG_FREE_ALL
#define LG_FREE_ALL LAGraph_Free ((void **) &mem, msg) ;

// w[index[i]] = min(w[index[i]], s[i]) for i in [0..n-1]
static GrB_Info Reduce_assign (GrB_Vector w,
        GrB_Vector s, GrB_Index *index, GrB_Index n, char *msg)
{
    GrB_Index *mem = NULL ;
    LG_TRY (LAGraph_Malloc ((void **) &mem, n*3, sizeof (GrB_Index), msg)) ;
    GrB_Index *ind = mem, *sval = mem + n, *wval = sval + n;
    LG_TRY (GrB_Vector_extractTuples(ind, wval, &n, w));
    LG_TRY (GrB_Vector_extractTuples(ind, sval, &n, s));
    for (GrB_Index i = 0; i < n; i++)
        if (sval[i] < wval[index[i]])
            wval[index[i]] = sval[i];
    LG_TRY (GrB_Vector_clear(w));
    LG_TRY (GrB_Vector_build(w, ind, wval, n, GrB_PLUS_UINT64));
    LG_FREE_ALL ;
    return GrB_SUCCESS;
}

//****************************************************************************
typedef struct
{
    uint64_t *data;     // array to malloc / free
    uint64_t *weight;   // minimum edge weight for each vertex
    uint64_t *parent;   // parent of each vertex in the spanning forest
    uint64_t *partner;  // partner vertex in the spanning forest
} MSF_context;
#define MSF_CONT    \
"typedef struct\n"              \
"{\n"                           \
"    uint64_t *data;     \n"   \
"    uint64_t *weight;   \n"   \
"    uint64_t *parent;   \n"   \
"    uint64_t *partner;  \n"   \
"} MSF_context;\n"
// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component

void selectEdge (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
{
    uint64_t *aij = (uint64_t*) x;
    (*z) = (thunk->weight[i] == *aij) && (thunk->parent[j] == thunk->partner[i]);
}
#define SELECTEDGE  \
"void selectEdge\n"                                                                 \
"(bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)\n"\
"{\n"                                                                               \
"    uint64_t *aij = (uint64_t*) x;\n"                                              \
"    (*z) = (thunk->weight[i] == *aij) && (thunk->parent[j] == thunk->partner[i]);\n"\
"}"

// edge removal:
// A(i, j) is removed when parent[i] == parent[j]

void removeEdge (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
{
    (*z) = (thunk->parent[i] != thunk->parent[j]);
}
#define REMOVEEDGE  \
"void removeEdge\n"                                                                         \
"(bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)\n"\
"{\n"                                                                               \
"    (*z) = (thunk->parent[i] != thunk->parent[j]);\n"                              \
"}"

//****************************************************************************

#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&S);                              \
    GrB_free (&T);                              \
    LAGraph_Free ((void **) &SI, msg);          \
    LAGraph_Free ((void **) &SJ, msg);          \
    LAGraph_Free ((void **) &SX, msg);          \
    LAGraph_Free ((void **) &context.data, msg);\
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
    GrB_free (&contx_type);             \
}

//****************************************************************************
int LAGraph_msf
(
    GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,       // input matrix
    bool sanitize,      // if true, ensure A is symmetric
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    LG_CLEAR_MSG ;
    MSF_context context = {NULL, NULL, NULL, NULL};
    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, i = NULL, t = NULL,
        edge = NULL, cedge = NULL, mask = NULL, index = NULL;
    GrB_Index *SI = NULL, *SJ = NULL, *SX = NULL;
    GrB_Type contx_type = NULL;
    GrB_BinaryOp comb = NULL;
    GrB_Semiring combMin = NULL;
    GrB_UnaryOp fst = NULL, snd = NULL;

    GrB_IndexUnaryOp s1 = NULL, s2 = NULL;
    if (result == NULL || A == NULL) return (GrB_NULL_POINTER) ;

    GrB_Index ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, A));
    GRB_TRY (GrB_Matrix_ncols (&ncols, A));
    if (n != ncols) return (GrB_DIMENSION_MISMATCH) ;

    if (sanitize)
    {
        // S = A+A'
        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n, n));
        GRB_TRY (GrB_eWiseAdd (S, 0, 0, GrB_PLUS_UINT64, A, A, GrB_DESC_T1));
    }
    else
    {
        // Use the input as-is, and assume it is GrB_UINT64 and symmetric
        GRB_TRY (GrB_Matrix_dup (&S, A));
    }

    GRB_TRY (GrB_Matrix_new (&T, GrB_UINT64, n, n));
    GRB_TRY (GrB_Vector_new (&t, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&f, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&i, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&edge, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&cedge, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n));
    GRB_TRY (GrB_Vector_new (&index, GrB_UINT64, n));

    // temporary arrays
    LG_TRY (LAGraph_Malloc ((void **) &SI, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SJ, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SX, 2*n, sizeof (GrB_Index), msg)) ;

    // global arrays
    LG_TRY (LAGraph_Malloc ((void **) &context.data, 3 * n, sizeof (uint64_t), msg)) ;
    context.parent  = context.data;
    context.partner = context.data + n;
    context.weight  = context.data + 2 * n;

    // prepare vectors
    for (GrB_Index i = 0; i < n; i++)
        context.parent[i] = i;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        f, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        f, NULL, NULL, GrB_ROWINDEX_INT64, f, (uint64_t) 0, NULL));
    GRB_TRY (GrB_assign (i, 0, 0, f, GrB_ALL, 0, 0));

    // semiring & monoid
    GrB_Index inf = ((uint64_t) INT_MAX << 32) ^ INT_MAX;
    GRB_TRY (GrB_BinaryOp_new (
        &comb, (GxB_binary_function) combine, 
        GrB_UINT64, GrB_UINT64, GrB_UINT64
    ));
    GRB_TRY (GrB_Semiring_new (&combMin, GrB_MIN_MONOID_UINT64, comb));
    GRB_TRY (GrB_UnaryOp_new (
        &fst, (GxB_unary_function) get_fst, GrB_UINT64, GrB_UINT64));
    GRB_TRY (GrB_UnaryOp_new (
        &snd, (GxB_unary_function) get_snd, GrB_UINT64, GrB_UINT64));

    // context type
    GRB_TRY (GxB_Type_new (
        &contx_type, sizeof (MSF_context), "MSF_context", MSF_CONT));
        
    // ops for GrB_select
    GRB_TRY(GxB_IndexUnaryOp_new (
        &s1, (GxB_index_unary_function) selectEdge, GrB_BOOL, GrB_UINT64, 
        contx_type, "selectEdge", SELECTEDGE
    ));
    GRB_TRY(GxB_IndexUnaryOp_new (
        &s2, (void *) removeEdge, GrB_BOOL, GrB_UINT64, contx_type, 
        "removeEdge", REMOVEEDGE
    ));

    // the main computation
    GrB_Index nvals, diff, ntuples = 0, num;
    GRB_TRY (GrB_Matrix_nvals (&nvals, S));
    for (int iters = 1; nvals > 0; iters++)
    {
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        GRB_TRY (GrB_assign (edge, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY (GrB_mxv (edge, 0, GrB_MIN_UINT64, combMin, S, f, 0));
        // cedge[u] = children's minimum edge  | if u is a root
        //          = (INT_MAX, u)             | otherwise
        GRB_TRY (GrB_assign (t, 0, 0, (uint64_t) INT_MAX, GrB_ALL, 0, 0));
        GRB_TRY (GrB_eWiseMult (cedge, 0, 0, comb, t, i, 0));
        LG_TRY (Reduce_assign (cedge, edge, context.parent, n, msg));
        // if (f[u] == u) f[u] := snd(cedge[u])  -- the index part of the edge
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, f, i, 0));
        GRB_TRY (GrB_apply (f, mask, GrB_SECOND_UINT64, snd, cedge, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL));
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, i, t, 0));
        GRB_TRY (GrB_assign (f, mask, GrB_MIN_UINT64, i, GrB_ALL, 0, 0));

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, i, f, 0));
        GRB_TRY (GrB_assign (cedge, mask, 0, inf, GrB_ALL, 0, 0));

        // 2. every vertex tries to know whether one of its edges is selected
        GRB_TRY (GrB_extract (t, 0, 0, cedge, context.parent, n, 0));
        GRB_TRY (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, edge, t, 0));

        // 3. each root picks a vertex from its children to generate the solution
        GRB_TRY (GrB_assign (index, 0, 0, n, GrB_ALL, 0, 0));
        GRB_TRY (GrB_assign (index, mask, 0, i, GrB_ALL, 0, 0));
        GRB_TRY (GrB_assign (t, 0, 0, n, GrB_ALL, 0, 0));
        LG_TRY (Reduce_assign (t, index, context.parent, n, msg));
        GRB_TRY (GrB_extract (index, 0, 0, t, context.parent, n, 0));
        GRB_TRY (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, i, index, 0));

        // 4. generate the select function (set the global pointers)
        GRB_TRY (GrB_assign (t, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY (GrB_apply (t, mask, 0, fst, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (NULL, context.weight, &n, t));
        GRB_TRY (GrB_assign (t, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY (GrB_apply (t, mask, 0, snd, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (NULL, context.partner, &n, t));
        GRB_TRY (GrB_Matrix_select_UDT (T, NULL, NULL, s1, S, &context, NULL));
        GRB_TRY (GrB_Vector_clear (t));

        // 5. the generated matrix may still have redundant edges
        //    remove the duplicates by GrB_mxv() and store them as tuples
        GRB_TRY (GrB_Vector_clear (edge));
        GRB_TRY (GrB_mxv (edge, mask, GrB_MIN_UINT64, combMin, T, i, 0));
        GRB_TRY (GrB_Vector_nvals (&num, edge));
        GRB_TRY (GrB_apply (t, 0, 0, snd, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (SI + ntuples, SJ + ntuples, &num, t));
        GRB_TRY (GrB_apply (t, 0, 0, fst, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (SI + ntuples, SX + ntuples, &num, t));
        GRB_TRY (GrB_Vector_clear (t));
        ntuples += num;

        // path halving until every vertex points on a root
        do {
            GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL));
            GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, f, t, 0));
            GrB_Vector temp = f;
            f = t;
            t = temp;
            temp = NULL;
            GRB_TRY (GrB_reduce (&diff, 0, GrB_PLUS_MONOID_UINT64, mask, 0));
        } while (diff != 0);

        // remove the edges in the same connected component
        GRB_TRY (GrB_Vector_extractTuples (NULL, context.parent, &n, f));
        GRB_TRY (GrB_Matrix_select_UDT (S, NULL, NULL, s2, S, &context, NULL)) ;
        GrB_Matrix_nvals (&nvals, S);
        if (nvals == 0) break;
    }
    GRB_TRY (GrB_Matrix_clear (T));
    GRB_TRY (GrB_Matrix_build (T, SI, SJ, SX, ntuples, GxB_IGNORE_DUP));
    *result = T;
    T = NULL ;

    LG_FREE_ALL;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
