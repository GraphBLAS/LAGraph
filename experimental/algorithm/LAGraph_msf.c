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

// TODO: a "sanitize" input is fine for now in the experimental folder, but it
// doesn't fit with the standard LAGraph API.  It will need to be removed when
// this method is moved to the src folder.  The input will also become an
// LAGraph_Graph, not a plain GrB_Matrix A.

#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

typedef struct
{
    union{
        int64_t wInt;
        double wFp;
    };
    uint64_t idx;
} pairW;

#define PAIRW               \
"typedef struct\n"          \
"{\n"                       \
"    union{\n"              \
"        int64_t wInt;\n"  \
"        double wFp;\n"     \
"    };\n"                  \
"    uint64_t idx;\n"       \
"} pairW;\n"

typedef struct
{
    uint64_t    *parent;   // parent of each vertex in the spanning forest
    pairW       *w_partner;  // partner vertex in the spanning forest
    // GrB_Type_Code type;
} MSF_context;

//pairW works inconsistently with JIT
#define MSF_CONT                \
"typedef struct\n"              \
"{\n"                           \
"    uint64_t   *parent;\n"     \
"    struct\n"                  \
"    {\n"                       \
"        union{\n"              \
"            int64_t wInt;\n"  \
"            double wFp;\n"     \
"        };\n"                  \
"        uint64_t idx;\n"       \
"    } *w_partner;\n"           \
"} MSF_context;\n"


//****************************************************************************
// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component

void selectEdge (bool *z, const int64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
{
    (*z) = (thunk->w_partner[i].wInt == *x) && (thunk->parent[j] == thunk->w_partner[i].idx);
}
#define SELECTEDGE  \
"void selectEdge\n"                                                                 \
"(bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)\n"\
"{\n"                                                                               \
"    (*z) = (thunk->w_partner[i].wInt == *x) && (thunk->parent[j] == thunk->w_partner[i].idx);\n"\
"}"

// edge removal:
// A(i, j) is removed when parent[i] == parent[j]

void removeEdge (bool *z, const int64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
{
    (*z) = (thunk->parent[i] != thunk->parent[j]);
}
#define REMOVEEDGE  \
"void removeEdge\n"                                                                 \
"(bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)\n"\
"{\n"                                                                               \
"    (*z) = (thunk->parent[i] != thunk->parent[j]);\n"                              \
"}"

//****************************************************************************

static void combine (pairW *z, const int64_t *x, const uint64_t *y)
{
    z->wInt = *x;
    z->idx = *y;
}
#define COMBINE \
"void combine (pairW *z, const uint64_t *x, const uint64_t *y)\n"\
"{\n"\
"    z->wInt = *x;\n"\
"    z->idx = *y;\n"\
"}\n"

static void get_fst (uint64_t *y, const pairW *x)
{
    *y = x->wInt;
}
#define GETFST \
"void get_fst (uint64_t *y, const pairW *x)\n"\
"{\n"\
"    *y = x->wInt;\n"\
"}\n"

static void get_snd (uint64_t *y, const pairW *x)
{
    *y = x->idx;
}
#define GETSND \
"void get_snd (uint64_t *y, const pairW *x)\n"\
"{\n"\
"    *y = x->idx;\n"\
"}\n"

static void tupleMinInt(pairW *z, const pairW *x, const pairW *y)
{
    bool xSmaller = x->wInt < y->wInt || (x->wInt == y->wInt && x->idx < y->idx);
    z->wInt = (xSmaller)? x->wInt: y->wInt;
    z->idx = (xSmaller)? x->idx: y->idx;
}
#define TUPLEMININT \
"void tupleMinInt(pairW *z, const pairW *x, const pairW *y)\n"\
"{\n"\
"    bool xSmaller = x->wInt < y->wInt || (x->wInt == y->wInt && x->idx < y->idx);\n"\
"    z->wInt = (xSmaller)? x->wInt: y->wInt; \n"\
"    z->idx = (xSmaller)? x->idx: y->idx; \n"\
"}\n"
static void tupleMinFp(pairW *z, const pairW *x, const pairW *y)
{
    bool xSmaller = x->wFp < y->wFp || (x->wFp == y->wFp && x->idx < y->idx);
    z->wFp = (xSmaller)? x->wFp: y->wFp;
    z->idx = (xSmaller)? x->idx: y->idx;
}
#define TUPLEMINFP \
"void tupleMinFp(pairW *z, const pairW *x, const pairW *y)\n"\
"{\n"\
"    bool xSmaller = x->wFp < y->wFp || (x->wFp == y->wFp && x->idx < y->idx);\n"\
"    z->wFp = (xSmaller)? x->wFp: y->wFp; \n"\
"    z->idx = (xSmaller)? x->idx: y->idx; \n"\
"}\n"

// Set z to the second -- sets bits regardless of weight type.
static void tuple2nd(pairW *z, const void *x, const pairW *y)
{
    z->wInt = y->wInt;
    z->idx = y->idx;
}
#define TUPLE2ND \
"void tuple2nd(pairW *z, const void *x, const pairW *y)\n"\
"{\n"\
"    z->wInt = y->wInt;\n"\
"    z->idx = y->idx;\n"\
"}\n"

// Since no arithmetic is done on tuples, can compare the ints to know if they
// are equal.
static void tupleEq(bool *z, const pairW *x, const pairW *y)
{
    *z = x->wInt == y->wInt && x->idx == y->idx;
}
#define TUPLEEQ \
"void tupleEq(bool *z, const pairW *x, const pairW *y)\n"\
"{\n"\
"    *z = x->wInt == y->wInt && x->idx == y->idx;\n"\
"}\n"

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                 \
{                                                   \
    GrB_free (&S);                                  \
    GrB_free (&T);                                  \
    LAGraph_Free ((void **) &SI, msg);              \
    LAGraph_Free ((void **) &SJ, msg);              \
    LAGraph_Free ((void **) &SX, msg);              \
    LAGraph_Free ((void **) &context.parent, msg);  \
    GrB_free (&f);                                  \
    GrB_free (&I);                                  \
    GrB_free (&t);                                  \
    GrB_free (&edge);                               \
    GrB_free (&cedge);                              \
    GrB_free (&tedge);                              \
    GrB_free (&mask);                               \
    GrB_free (&index_v);                            \
    GrB_free (&comb);                               \
    GrB_free (&minComb);                            \
    GrB_free (&fst);                                \
    GrB_free (&snd);                                \
    GrB_free (&s1);                                 \
    GrB_free (&s2);                                 \
    GrB_free (&contx_type);                         \
    GrB_free (&parent_v);                           \
    GrB_free (&ramp);                               \
    GrB_free (&pairMin);                            \
    GrB_free (&pairSec);                            \
    GrB_free (&pairEq);                             \
    GrB_free (&pairMin_monoid);                     \
    GrB_free (&pairMin2nd);                         \
    GrB_free (&lg_pair);                            \
    GrB_free (&max_weight);                         \
}
int LAGraph_msf
(
    GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,       // input matrix
    bool sanitize,      // if true, ensure A is symmetric
    char *msg
)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LG_CLEAR_MSG ;
    MSF_context context = {
        .parent = NULL, .w_partner = NULL,
        // .type = GrB_UINT64_CODE
    };
    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, I = NULL, t = NULL, parent_v = NULL, tedge = NULL,
        edge = NULL, cedge = NULL, mask = NULL, index_v = NULL, ramp = NULL;

    GrB_Index *SI = NULL, *SJ = NULL;
    void *SX = NULL;
    GrB_Type contx_type = NULL, lg_pair = NULL, weight_type = NULL;
    GrB_BinaryOp comb = NULL, pairMin = NULL, pairSec = NULL, pairEq = NULL;
    GrB_Monoid pairMin_monoid = NULL;
    GrB_Semiring minComb = NULL, pairMin2nd = NULL;
    GrB_UnaryOp fst = NULL, snd = NULL;
    GrB_Scalar max_weight = NULL;
    int edge_h = GrB_DEFAULT;
    uint64_t edge_size = 0, edge_n = 0;
    GrB_IndexUnaryOp s1 = NULL, s2 = NULL;


    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    if (result == NULL || A == NULL) return (GrB_NULL_POINTER) ;
    GrB_Index ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    LG_ASSERT(n == ncols, GrB_DIMENSION_MISMATCH) ;

    GrB_Type_Code tcode = 0;
    if (sanitize)
    {
        GrB_Matrix_get_INT32(A, (int *) &(tcode), GrB_EL_TYPE_CODE);
        // S = A+A'
        switch (tcode)
        {
            case GrB_INT8_CODE:
            case GrB_INT16_CODE:
            case GrB_INT32_CODE:
            case GrB_INT64_CODE:
            case GrB_BOOL_CODE:
            case GrB_UINT8_CODE:
            case GrB_UINT16_CODE:
            case GrB_UINT32_CODE:
            case GrB_UINT64_CODE:
                tcode = GrB_INT64_CODE;
                GRB_TRY (GrB_Matrix_new (&S, GrB_INT64, n, n)) ;
                GRB_TRY (GrB_Matrix_eWiseAdd_BinaryOp
                    (S, NULL, NULL, GrB_MIN_INT64, A, A, GrB_DESC_T1)) ;
                break;
            case GrB_FP32_CODE:
            case GrB_FP64_CODE:
                tcode = GrB_FP64_CODE;
                GRB_TRY (GrB_Matrix_new (&S, GrB_FP64, n, n)) ;
                GRB_TRY (GrB_Matrix_eWiseAdd_BinaryOp
                    (S, NULL, NULL, GrB_MIN_FP64, A, A, GrB_DESC_T1)) ;
                break;
            default:
                LG_ASSERT(false, GrB_DOMAIN_MISMATCH) ;
                break;
        }
        weight_type = (tcode == GrB_INT64_CODE) ? GrB_INT64 : GrB_FP64 ;
    }
    else
    {
        // Use the input as-is, and assume it is symmetric
        GrB_Matrix_get_INT32(A, (int *) &(tcode), GrB_EL_TYPE_CODE) ;
        LG_ASSERT(tcode < 12 && tcode > 0, GrB_DOMAIN_MISMATCH) ;
        tcode = (tcode == GrB_FP32_CODE || tcode == GrB_FP64_CODE)?
            GrB_FP64_CODE: GrB_INT64_CODE;
        weight_type = (tcode == GrB_INT64_CODE) ? GrB_INT64 : GrB_FP64 ;
        GRB_TRY (GrB_Matrix_new (&S, weight_type, n, n)) ;
        GRB_TRY (GrB_Matrix_assign
                (S, NULL, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL)) ;
    }
    GRB_TRY (GxB_Type_new   (&lg_pair, sizeof(pairW), "pairW", PAIRW)) ;
    GRB_TRY (GrB_Matrix_new (&T, weight_type, n, n)) ;
    GRB_TRY (GrB_Vector_new (&t, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&f, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&ramp, GrB_INT64, n + 1)) ;
    GRB_TRY (GrB_Vector_new (&edge, lg_pair, n)) ;
    GRB_TRY (GrB_Vector_new (&cedge, lg_pair, n)) ;
    GRB_TRY (GrB_Vector_new (&tedge, lg_pair, n)) ;
    GRB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&index_v, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&parent_v, GrB_UINT64, n)) ;

    LG_TRY (LAGraph_Malloc  ((void **) &SI, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc  ((void **) &SJ, 2*n, sizeof (GrB_Index), msg)) ;
    size_t sx_size = (tcode == GrB_INT64_CODE) ? sizeof (int64_t) : sizeof (double) ;
    LG_TRY (LAGraph_Malloc  (&SX, 2*n, sx_size, msg)) ;

    // context arrays
    LG_TRY (LAGraph_Malloc
        ((void **) &context.parent, n, sizeof (uint64_t), msg)) ;

    // prepare vectors
    for (uint64_t i = 0; i < n; i++)
        context.parent[i] = i;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        f, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        f, NULL, NULL, GrB_ROWINDEX_INT64, f, (int64_t) 0, NULL)) ;
    GRB_TRY (GrB_Vector_dup (&I, f)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        ramp, NULL, NULL, (uint64_t) 0, GrB_ALL, n + 1, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        ramp, NULL, NULL, GrB_ROWINDEX_INT64, ramp, (int64_t) 0, NULL)) ;
    GRB_TRY (GxB_Vector_load(parent_v, (void **) &context.parent,
        GrB_UINT64, n, 3 * n * sizeof (uint64_t), GxB_IS_READONLY, NULL)) ;
    // semiring & monoid
    pairW inf = {.wInt = INT64_MAX, .idx = UINT64_MAX};
    if(tcode == GrB_FP64_CODE) inf.wFp = INFINITY;

    GRB_TRY (GxB_Scalar_new(&max_weight, weight_type)) ;


    GRB_TRY (GxB_BinaryOp_new (
        &comb, (GxB_binary_function) combine,
        lg_pair, weight_type, GrB_UINT64, "combine", COMBINE
    )) ;

    if(tcode == GrB_INT64_CODE)
    {
        GRB_TRY (GxB_Scalar_setElement_INT64(max_weight, INT64_MAX)) ;
        GRB_TRY (GxB_BinaryOp_new (
            &pairMin, (GxB_binary_function) tupleMinInt,
            lg_pair, lg_pair, lg_pair, "tupleMinInt", TUPLEMININT
        )) ;
    }
    else
    {
        GRB_TRY (GxB_Scalar_setElement_FP64(max_weight, INFINITY)) ;
        GRB_TRY (GxB_BinaryOp_new (
            &pairMin, (GxB_binary_function) tupleMinFp,
            lg_pair, lg_pair, lg_pair, "tupleMinFp", TUPLEMINFP
        )) ;
    }

    GRB_TRY (GxB_BinaryOp_new (
        &pairSec, (GxB_binary_function) tuple2nd,
        lg_pair, GrB_BOOL, lg_pair, "tuple2nd", TUPLE2ND
    )) ;

    GRB_TRY (GxB_BinaryOp_new (
        &pairEq, (GxB_binary_function) tupleEq,
        GrB_BOOL, lg_pair, lg_pair, "tupleEq", TUPLEEQ
    )) ;
    
    GRB_TRY (GrB_Monoid_new_UDT (&pairMin_monoid, pairMin, (void *) &inf)) ;
    GRB_TRY (GrB_Semiring_new (&minComb, pairMin_monoid, comb)) ;
    GRB_TRY (GrB_Semiring_new (&pairMin2nd, pairMin_monoid, pairSec)) ;

    GRB_TRY (GxB_UnaryOp_new (
        &fst, (GxB_unary_function) get_fst, weight_type, lg_pair,
        "get_fst", GETFST)) ;
    GRB_TRY (GxB_UnaryOp_new (
        &snd, (GxB_unary_function) get_snd, GrB_UINT64, lg_pair,
        "get_snd", GETSND)) ;

    // context type
    GRB_TRY (GxB_Type_new (
        &contx_type, sizeof (MSF_context), "MSF_context", MSF_CONT)) ;

    // ops for GrB_select
    GRB_TRY(GxB_IndexUnaryOp_new (
        &s1, (GxB_index_unary_function) selectEdge, GrB_BOOL, weight_type,
        contx_type, "selectEdge", SELECTEDGE
    )) ;
    GRB_TRY(GxB_IndexUnaryOp_new (
        &s2, (void *) removeEdge, GrB_BOOL, GrB_UINT64, contx_type,
        "removeEdge", REMOVEEDGE
    )) ;

    // the main computation
    GrB_Index nvals, ntuples = 0, num;
    bool diff = false;
    GRB_TRY (GrB_Matrix_nvals (&nvals, S)) ;
    for (int iters = 1; nvals > 0; iters++)
    {
        #ifdef DEBUG
        LG_ASSERT(iters < 100, LAGRAPH_CONVERGENCE_FAILURE);
        #endif
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        GRB_TRY (GrB_Vector_assign_UDT (
            edge, NULL, NULL, (void *) &inf, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_mxv (edge, NULL, pairMin, minComb, S, f, NULL)) ;

        // cedge[u] = children's minimum edge  | if u is a root
        //          = (max_weight, u)          | otherwise
        GRB_TRY (GrB_Vector_apply_BinaryOp1st_Scalar (
            cedge, NULL, NULL, comb, max_weight, I, NULL)) ;
        LG_TRY (LAGraph_FastAssign_Semiring(
            cedge, NULL, pairMin, parent_v, edge, ramp, pairMin2nd, NULL, msg
        )) ;
        // if (f[u] == u) f[u] := snd(cedge[u])  -- the index part of the edge
        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, f, I, NULL)) ;
        GRB_TRY (GrB_apply (f, mask, GrB_SECOND_UINT64, snd, cedge, NULL)) ;
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL)) ;
        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, I, t, NULL)) ;
        GRB_TRY (GrB_assign (f, mask, GrB_MIN_UINT64, I, GrB_ALL, 0, NULL)) ;

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, I, f, NULL)) ;
        GRB_TRY (GrB_assign (cedge, mask, NULL, (void *) &inf, GrB_ALL, 0, NULL)) ;

        // 2. every vertex tries to know whether one of its edges is selected
        GRB_TRY (GxB_Vector_extract_Vector (
            tedge, NULL, NULL, cedge, parent_v, NULL)) ;
        GRB_TRY (GrB_eWiseMult (mask ,NULL, NULL, pairEq, edge, tedge, NULL)) ;

        // 3. each root picks a vertex from its children to generate the solution
        GRB_TRY (GrB_assign (index_v, NULL, NULL, n, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_assign (index_v, mask, NULL, I, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_assign (t, NULL, NULL, n, GrB_ALL, 0, NULL)) ;
        LG_TRY (LAGraph_FastAssign_Semiring(
            t, NULL, GrB_MIN_UINT64, parent_v, index_v, ramp,
            GrB_MIN_SECOND_SEMIRING_UINT64, NULL, msg
        )) ;
        GRB_TRY (GxB_Vector_extract_Vector (
            index_v, NULL, NULL, t, parent_v, NULL)) ;
        GRB_TRY (GrB_eWiseMult (mask ,NULL, NULL, GrB_EQ_UINT64, I, index_v, NULL)) ;

        // 4. generate the select function (set the global pointers)
        GRB_TRY (GxB_Vector_unload(
            edge, (void **) &context.w_partner, &lg_pair, &edge_n, &edge_size,
            &edge_h, NULL)) ;
        GRB_TRY (GrB_Matrix_select_UDT (T, NULL, NULL, s1, S, &context, NULL)) ;
        GRB_TRY (GxB_Vector_load(
            edge, (void **) &context.w_partner, lg_pair, edge_n, edge_size,
            edge_h, NULL)) ;
        GRB_TRY (GrB_Vector_clear (t)) ;

        // 5. the generated matrix may still have redundant edges
        //    remove the duplicates by GrB_mxv() and store them as tuples
        GRB_TRY (GrB_Vector_clear (edge)) ;
        GRB_TRY (GrB_mxv (edge, mask, pairMin, minComb, T, I, NULL)) ;
        GRB_TRY (GrB_Vector_nvals (&num, edge)) ;
        GRB_TRY (GrB_apply (t, NULL, NULL, snd, edge, NULL)) ;
        GRB_TRY (GrB_Vector_extractTuples (NULL, SJ + ntuples, &num, t)) ;
        if(tcode == GrB_INT64_CODE)
        {
            GRB_TRY (GrB_apply (t, NULL, NULL, fst, edge, NULL)) ;
            GRB_TRY (GrB_Vector_extractTuples_INT64 (
                SI + ntuples, ((int64_t *) SX) + ntuples, &num, t)) ;
            GRB_TRY (GrB_Vector_clear (t)) ;
        }
        else
        {
            GRB_TRY (GrB_free(&t)) ;
            GRB_TRY (GrB_Vector_new(&t, weight_type, n)) ;
            GRB_TRY (GrB_apply (t, NULL, NULL, fst, edge, NULL)) ;
            GRB_TRY (GrB_Vector_extractTuples_FP64 (
                SI + ntuples, ((double *) SX) + ntuples, &num, t)) ;
            GRB_TRY (GrB_Vector_clear (t)) ;
            GRB_TRY (GrB_free(&t)) ;
            GRB_TRY (GrB_Vector_new(&t, GrB_UINT64, n)) ;
        }

        ntuples += num;

        // path halving until every vertex points on a root
        do {
            GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL)) ;
            GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_NE_UINT64, f, t, NULL)) ;
            GrB_Vector temp = f;
            f = t;
            t = temp;
            GRB_TRY (GrB_Vector_reduce_BOOL (&diff, NULL, GrB_LOR_MONOID_BOOL, mask, NULL)) ;
        } while (diff);

        // remove the edges in the same connected component
        GRB_TRY (GrB_Vector_extractTuples (NULL, context.parent, &n, f)) ;
        GRB_TRY (GrB_Matrix_select_UDT (S, NULL, NULL, s2, S, &context, NULL)) ;
        GrB_Matrix_nvals (&nvals, S);
        if (nvals == 0) break;
    }

    GRB_TRY (GrB_Matrix_clear (T)) ;
    if(tcode == GrB_INT64_CODE)
    {
        GRB_TRY (GrB_Matrix_build_INT64 (
            T, SI, SJ, (int64_t *)SX, ntuples, GxB_IGNORE_DUP)) ;
    }
    else
    {
        GRB_TRY (GrB_Matrix_build_FP64 (
            T, SI, SJ, (double *)SX, ntuples, GxB_IGNORE_DUP)) ;
    }

    *result = T;
    T = NULL ;

    LG_FREE_ALL;
    return (GrB_SUCCESS) ;
    #else
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}
