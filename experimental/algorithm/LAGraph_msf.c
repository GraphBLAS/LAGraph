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

#if LG_SUITESPARSE_GRAPHBLAS_V10
typedef struct
{
    union{
        uint64_t wInt;
        double wFp;
    };
    uint64_t idx;
} pairW;

#define PAIRW               \
"typedef struct\n"          \
"{\n"                       \
"    union{\n"              \
"        uint64_t wInt;\n"  \
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
"            uint64_t wInt;\n"  \
"            double wFp;\n"     \
"        };\n"                  \
"        uint64_t idx;\n"       \
"    } *w_partner;\n"           \
/* "    GrB_Type_Code type;\n"          */\
"} MSF_context;\n"      


//****************************************************************************
// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component

void selectEdge (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
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

static void combine (pairW *z, const uint64_t *x, const uint64_t *y)
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
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&S);                              \
    GrB_free (&T);                              \
    LAGraph_Free ((void **) &SI, msg);          \
    LAGraph_Free ((void **) &SJ, msg);          \
    LAGraph_Free ((void **) &SX, msg);          \
    LAGraph_Free ((void **) &context.parent, msg);\
    GrB_free (&f);                      \
    GrB_free (&I);                      \
    GrB_free (&t);                      \
    GrB_free (&edge);                   \
    GrB_free (&cedge);                  \
    GrB_free (&tedge);                  \
    GrB_free (&mask);                   \
    GrB_free (&index);                  \
    GrB_free (&comb);                   \
    GrB_free (&minComb);                \
    GrB_free (&fst);                    \
    GrB_free (&snd);                    \
    GrB_free (&s1);                     \
    GrB_free (&s2);                     \
    GrB_free (&contx_type);             \
    GrB_free (&parent_v);               \
    GrB_free (&ramp);                   \
    GrB_free (&pairMin);                \
    GrB_free (&pairSec);                \
    GrB_free (&pairEq);                 \
    GrB_free (&pairMin_monoid);         \
    GrB_free (&pairMin2nd);             \
    GrB_free (&lg_pair);                \
}
int LAGraph_msf_GB10
(
    GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,       // input matrix
    bool sanitize,      // if true, ensure A is symmetric
    char *msg
)
{
    LG_CLEAR_MSG ;
    MSF_context context = {
        .parent = NULL, .w_partner = NULL, 
        // .type = GrB_UINT64_CODE
    };
    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, I = NULL, t = NULL, parent_v = NULL, tedge = NULL,
        edge = NULL, cedge = NULL, mask = NULL, index = NULL, ramp = NULL;

    GrB_Index *SI = NULL, *SJ = NULL, *SX = NULL;
    GrB_Type contx_type = NULL, lg_pair = NULL, weight_type = NULL;
    GrB_BinaryOp comb = NULL, pairMin = NULL, pairSec = NULL, pairEq = NULL;
    GrB_Monoid pairMin_monoid = NULL;
    GrB_Semiring minComb = NULL, pairMin2nd = NULL;
    GrB_UnaryOp fst = NULL, snd = NULL;
    int edge_h = GrB_DEFAULT;
    uint64_t edge_size = 0, edge_n = 0; 
    GrB_IndexUnaryOp s1 = NULL, s2 = NULL;


    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    
    if (result == NULL || A == NULL) return (GrB_NULL_POINTER) ;
    GrB_Index ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, A));
    GRB_TRY (GrB_Matrix_ncols (&ncols, A));
    LG_ASSERT(n == ncols, GrB_DIMENSION_MISMATCH) ;

    GrB_Type_Code tcode = 0;
    if (sanitize)
    {
        GrB_Matrix_get_INT32(A, &(tcode), GrB_EL_TYPE_CODE);
        // S = A+A'
        switch (tcode)
        {
        case GrB_INT8_CODE:
        case GrB_INT16_CODE:
        case GrB_INT32_CODE:
        case GrB_INT64_CODE:
            // TODO handle negative weights here.
        case GrB_BOOL_CODE:
        case GrB_UINT8_CODE:
        case GrB_UINT16_CODE:
        case GrB_UINT32_CODE:
        case GrB_UINT64_CODE:
            tcode = GrB_UINT64_CODE;
            GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n, n));
            GRB_TRY (GrB_Matrix_eWiseAdd_BinaryOp 
                (S, NULL, NULL, GrB_MIN_UINT64, A, A, GrB_DESC_T1));
            break;
        case GrB_FP32_CODE:
        case GrB_FP64_CODE:
            tcode = GrB_FP64_CODE;
            GRB_TRY (GrB_Matrix_new (&S, GrB_FP64, n, n));
            GRB_TRY (GrB_Matrix_eWiseAdd_BinaryOp 
                (S, NULL, NULL, GrB_MIN_FP64, A, A, GrB_DESC_T1));
            break;
        default:
            LG_ASSERT(false, GrB_DOMAIN_MISMATCH) ;
            break;
        }
        weight_type = tcode == GrB_UINT64_CODE? GrB_UINT64: GrB_FP64;
    }
    else
    {
        // Use the input as-is, and assume it is symmetric
        GrB_Matrix_get_INT32(A, &(tcode), GrB_EL_TYPE_CODE);
        LG_ASSERT(tcode < 12 && tcode > 0, GrB_DOMAIN_MISMATCH);
        tcode = (tcode == GrB_FP32_CODE || tcode == GrB_FP64_CODE)? 
            GrB_FP64_CODE: GrB_UINT64_CODE;
        weight_type = tcode == GrB_UINT64_CODE? GrB_UINT64: GrB_FP64;
        GRB_TRY (GrB_Matrix_new (&S, weight_type, n, n));
        GRB_TRY (GrB_Matrix_assign
                (S, NULL, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL));
    }
    GRB_TRY (GxB_Type_new(&lg_pair, sizeof(pairW), "pairW", PAIRW));
    GRB_TRY (GrB_Matrix_new (&T, weight_type, n, n));
    GRB_TRY (GrB_Vector_new (&t, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&f, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&ramp, GrB_INT64, n + 1));
    GRB_TRY (GrB_Vector_new (&edge, lg_pair, n));
    GRB_TRY (GrB_Vector_new (&cedge, lg_pair, n));
    GRB_TRY (GrB_Vector_new (&tedge, lg_pair, n));
    GRB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n));
    GRB_TRY (GrB_Vector_new (&index, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&parent_v, GrB_UINT64, n));

    LG_TRY (LAGraph_Malloc ((void **) &SI, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SJ, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SX, 2*n, sizeof (GrB_Index), msg)) ;

    // context arrays
    LG_TRY (LAGraph_Malloc ((void **) &context.parent, n, sizeof (uint64_t), msg)) ;

    // prepare vectors
    for (uint64_t i = 0; i < n; i++)
        context.parent[i] = i;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        f, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        f, NULL, NULL, GrB_ROWINDEX_INT64, f, (int64_t) 0, NULL));
    GRB_TRY (GrB_Vector_dup (&I, f));
    GRB_TRY (GrB_Vector_assign_UINT64 (
        ramp, NULL, NULL, (uint64_t) 0, GrB_ALL, n + 1, NULL));
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        ramp, NULL, NULL, GrB_ROWINDEX_INT64, ramp, (int64_t) 0, NULL));
    GRB_TRY (GxB_Vector_load(parent_v, (void **) &context.parent, 
        GrB_UINT64, n, 3 * n * sizeof (uint64_t), GxB_IS_READONLY, NULL));
    // semiring & monoid
    pairW inf = {.wInt = UINT64_MAX, .idx = UINT64_MAX};
    if(tcode == GrB_FP64_CODE) inf.wFp = INFINITY;

    GRB_TRY (GxB_BinaryOp_new (
        &comb, (GxB_binary_function) combine,
        lg_pair, weight_type, GrB_UINT64, "combine", COMBINE
    ));
    if(tcode == GrB_UINT64_CODE)
    {
        GRB_TRY (GxB_BinaryOp_new (
            &pairMin, (GxB_binary_function) tupleMinInt, 
            lg_pair, lg_pair, lg_pair, "tupleMinInt", TUPLEMININT
        ));
    }
    else
    {
        GRB_TRY (GxB_BinaryOp_new (
            &pairMin, (GxB_binary_function) tupleMinFp, 
            lg_pair, lg_pair, lg_pair, "tupleMinFp", TUPLEMINFP
        ));
    }
    
    GRB_TRY (GxB_BinaryOp_new (
        &pairSec, (GxB_binary_function) tuple2nd, 
        lg_pair, GrB_BOOL, lg_pair, "tuple2nd", TUPLE2ND
    ));
    GRB_TRY (GxB_BinaryOp_new (
        &pairEq, (GxB_binary_function) tupleEq, 
        GrB_BOOL, lg_pair, lg_pair, "tupleEq", TUPLEEQ
    ));
    GRB_TRY (GrB_Monoid_new_UDT (&pairMin_monoid, pairMin, (void *) &inf));
    GRB_TRY (GrB_Semiring_new (&minComb, pairMin_monoid, comb));
    GRB_TRY (GrB_Semiring_new (&pairMin2nd, pairMin_monoid, pairSec));
    GRB_TRY (GxB_UnaryOp_new (
        &fst, (GxB_unary_function) get_fst, weight_type, lg_pair, 
        "get_fst", GETFST));
    GRB_TRY (GxB_UnaryOp_new (
        &snd, (GxB_unary_function) get_snd, GrB_UINT64, lg_pair,
        "get_snd", GETSND));

    // context type
    GRB_TRY (GxB_Type_new (
        &contx_type, sizeof (MSF_context), "MSF_context", MSF_CONT));
        
    // ops for GrB_select
    GRB_TRY(GxB_IndexUnaryOp_new (
        &s1, (GxB_index_unary_function) selectEdge, GrB_BOOL, weight_type, 
        contx_type, "selectEdge", SELECTEDGE
    ));
    GRB_TRY(GxB_IndexUnaryOp_new (
        &s2, (void *) removeEdge, GrB_BOOL, GrB_UINT64, contx_type, 
        "removeEdge", REMOVEEDGE
    ));

    // the main computation
    GrB_Index nvals, ntuples = 0, num;
    bool diff = false;
    GRB_TRY (GrB_Matrix_nvals (&nvals, S));
    for (int iters = 1; nvals > 0; iters++)
    {
        #ifdef DEBUG
        LG_ASSERT(iters < 100, LAGRAPH_CONVERGENCE_FAILURE);
        #endif
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        GRB_TRY (GrB_Vector_assign_UDT (
            edge, NULL, NULL, (void *) &inf, GrB_ALL, 0, NULL));
        GRB_TRY (GrB_mxv (edge, 0, pairMin, minComb, S, f, 0));
        // cedge[u] = children's minimum edge  | if u is a root
        //          = (INT_MAX, u)             | otherwise
        GRB_TRY (GrB_assign (t, 0, 0, (uint64_t) INT_MAX, GrB_ALL, 0, 0));
        GRB_TRY (GrB_eWiseMult (cedge, 0, 0, comb, t, I, 0));
        LG_TRY (LAGraph_FastAssign_Semiring(
            cedge, NULL, pairMin, parent_v, edge, ramp, pairMin2nd, NULL, msg
        ));
        // if (f[u] == u) f[u] := snd(cedge[u])  -- the index part of the edge
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, f, I, 0));
        GRB_TRY (GrB_apply (f, mask, GrB_SECOND_UINT64, snd, cedge, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL));
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, I, t, 0));
        GRB_TRY (GrB_assign (f, mask, GrB_MIN_UINT64, I, GrB_ALL, 0, 0));

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, I, f, 0));
        GRB_TRY (GrB_assign (cedge, mask, 0, (void *) &inf, GrB_ALL, 0, 0));

        // 2. every vertex tries to know whether one of its edges is selected
        GRB_TRY (GxB_Vector_extract_Vector (
            tedge, NULL, NULL, cedge, parent_v, NULL));
        GRB_TRY (GrB_eWiseMult (mask ,0, 0, pairEq, edge, tedge, 0));

        // 3. each root picks a vertex from its children to generate the solution
        GRB_TRY (GrB_assign (index, 0, 0, n, GrB_ALL, 0, 0));
        GRB_TRY (GrB_assign (index, mask, 0, I, GrB_ALL, 0, 0));
        GRB_TRY (GrB_assign (t, 0, 0, n, GrB_ALL, 0, 0));
        LG_TRY (LAGraph_FastAssign_Semiring(
            t, NULL, GrB_MIN_UINT64, parent_v, index, ramp, 
            GrB_MIN_SECOND_SEMIRING_UINT64, NULL, msg
        ));
        GRB_TRY (GxB_Vector_extract_Vector (
            index, NULL, NULL, t, parent_v, NULL));
        GRB_TRY (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, I, index, 0));

        // 4. generate the select function (set the global pointers)
        GRB_TRY (GxB_Vector_unload(
            edge, (void **) &context.w_partner, &lg_pair, &edge_n, &edge_size,
            &edge_h, NULL)) ;
        GRB_TRY (GrB_Matrix_select_UDT (T, NULL, NULL, s1, S, &context, NULL));
        GRB_TRY (GxB_Vector_load(
            edge, (void **) &context.w_partner, lg_pair, edge_n, edge_size,
            edge_h, NULL)) ;
        GRB_TRY (GrB_Vector_clear (t));

        // 5. the generated matrix may still have redundant edges
        //    remove the duplicates by GrB_mxv() and store them as tuples
        GRB_TRY (GrB_Vector_clear (edge));
        GRB_TRY (GrB_mxv (edge, mask, pairMin, minComb, T, I, 0));
        GRB_TRY (GrB_Vector_nvals (&num, edge));
        GRB_TRY (GrB_apply (t, 0, 0, snd, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (NULL, SJ + ntuples, &num, t));
        if(tcode == GrB_UINT64_CODE)
        {
            GRB_TRY (GrB_apply (t, 0, 0, fst, edge, 0));
            GRB_TRY (GrB_Vector_extractTuples_UINT64 (
                SI + ntuples, SX + ntuples, &num, t));
            GRB_TRY (GrB_Vector_clear (t));
        }
        else
        {
            GRB_TRY (GrB_free(&t)) ;
            GRB_TRY (GrB_Vector_new(&t, weight_type, n));
            GRB_TRY (GrB_apply (t, 0, 0, fst, edge, 0));
            GRB_TRY (GrB_Vector_extractTuples_FP64 (
                SI + ntuples, (double *) SX + ntuples, &num, t));
            GRB_TRY (GrB_Vector_clear (t));
            GRB_TRY (GrB_free(&t)) ;
            GRB_TRY (GrB_Vector_new(&t, GrB_UINT64, n));
        }
        
        ntuples += num;

        // path halving until every vertex points on a root
        do {
            GRB_TRY (GxB_Vector_extract_Vector (t, NULL, NULL, f, f, NULL));
            GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, f, t, 0));
            GrB_Vector temp = f;
            f = t;
            t = temp;
            GRB_TRY (GrB_Vector_reduce_BOOL (&diff, NULL, GrB_LOR_MONOID_BOOL, mask, 0));
        } while (diff);

        // remove the edges in the same connected component
        GRB_TRY (GrB_Vector_extractTuples (NULL, context.parent, &n, f));
        GRB_TRY (GrB_Matrix_select_UDT (S, NULL, NULL, s2, S, &context, NULL)) ;
        GrB_Matrix_nvals (&nvals, S);
        if (nvals == 0) break;
    }
    GRB_TRY (GrB_Matrix_clear (T));
    if(tcode == GrB_UINT64_CODE)
    {
        GRB_TRY (GrB_Matrix_build_UINT64 (T, SI, SJ, (uint64_t *)SX, ntuples, GxB_IGNORE_DUP));
    }
    else
    {
        GRB_TRY (GrB_Matrix_build_FP64 (T, SI, SJ, (double *)SX, ntuples, GxB_IGNORE_DUP));
    }
    *result = T;
    T = NULL ;

    LG_FREE_ALL;
    return (GrB_SUCCESS) ;
}

#else
//****************************************************************************
// encode each edge into a single uint64_t
static void combine (uint64_t *z, const uint64_t *x, const uint64_t *y)
{
    *z = ((*x) << 32) | (*y);
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

typedef struct
{
    uint64_t *data;     // array to malloc / free
    uint64_t *parent;   // parent of each vertex in the spanning forest
    uint64_t *partner;  // partner vertex in the spanning forest
    uint64_t *weight;
    GrB_Type_Code type;
} MSF_context;
#define MSF_CONT    \
"typedef struct\n"              \
"{\n"                           \
"    uint64_t *data;     \n"   \
"    uint64_t *parent;   \n"   \
"    uint64_t *partner;  \n"   \
"    uint64_t *weight;   \n"   \
"} MSF_context;\n"

//****************************************************************************
// generate solution:
// for each element A(i, j), it is selected if
//   1. weight[i] == A(i, j)    -- where weight[i] stores i's minimum edge weight
//   2. parent[j] == partner[i] -- j belongs to the specified connected component

void selectEdge (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
{
    (*z) = (thunk->weight[i] == *x) && (thunk->parent[j] == thunk->partner[i]);
}
#define SELECTEDGE  \
"void selectEdge\n"                                                                 \
"(bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)\n"\
"{\n"                                                                               \
"    (*z) = (thunk->weight[i] == *x) && (thunk->parent[j] == thunk->partner[i]);\n"\
"}"

// edge removal:
// A(i, j) is removed when parent[i] == parent[j]

void removeEdge (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const MSF_context *thunk)
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

#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&S);                              \
    GrB_free (&T);                              \
    LAGraph_Free ((void **) &V, msg);           \
    LAGraph_Free ((void **) &SI, msg);          \
    LAGraph_Free ((void **) &SJ, msg);          \
    LAGraph_Free ((void **) &SX, msg);          \
    LAGraph_Free ((void **) &context.data, msg);\
    GrB_free (&f);                      \
    GrB_free (&I);                      \
    GrB_free (&t);                      \
    GrB_free (&edge);                   \
    GrB_free (&cedge);                  \
    GrB_free (&mask);                   \
    GrB_free (&index);                  \
    GrB_free (&comb);                   \
    GrB_free (&minComb);                \
    GrB_free (&fst);                    \
    GrB_free (&snd);                    \
    GrB_free (&s1);                     \
    GrB_free (&s2);                     \
    GrB_free (&contx_type);             \
    GrB_free (&parent_v);               \
    GrB_free (&ramp);                   \
}
#endif

//****************************************************************************
int LAGraph_msf
(
    GrB_Matrix *result, // output: an unsymmetrical matrix, the spanning forest
    GrB_Matrix A,       // input matrix
    bool sanitize,      // if true, ensure A is symmetric
    char *msg
)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10
    return LAGraph_msf_GB10(result, A, sanitize, msg);
#elif LAGRAPH_SUITESPARSE
    LG_CLEAR_MSG ;
    MSF_context context = {NULL, NULL, NULL, NULL};
    GrB_Info info;
    GrB_Index n;
    GrB_Matrix S = NULL, T = NULL;
    GrB_Vector f = NULL, I = NULL, t = NULL, parent_v = NULL,
        edge = NULL, cedge = NULL, mask = NULL, index = NULL, ramp = NULL;

    GrB_Index *SI = NULL, *SJ = NULL, *SX = NULL, *V = NULL;
    GrB_Type contx_type = NULL;
    GrB_BinaryOp comb = NULL;
    GrB_Semiring minComb = NULL;
    GrB_UnaryOp fst = NULL, snd = NULL;
    GrB_Type_Code tcode = 0;
    GrB_IndexUnaryOp s1 = NULL, s2 = NULL;
    if (result == NULL || A == NULL) return (GrB_NULL_POINTER) ;

    GrB_Index ncols ;
    GRB_TRY (GrB_Matrix_nrows (&n, A));
    GRB_TRY (GrB_Matrix_ncols (&ncols, A));
    LG_ASSERT(n == ncols, GrB_DIMENSION_MISMATCH) ;

    // // TODO: use tuple types so that we don't overflow.
    LG_ASSERT_MSG (n <= INT32_MAX, 
        GrB_INVALID_VALUE, "Matrix too large. Exiting to prevent overflow") ;
    uint64_t max_val = 0;
    GRB_TRY(GrB_Matrix_reduce_UINT64(
        &max_val, NULL, GrB_MAX_MONOID_UINT64, A, NULL)) ;
    LG_ASSERT_MSG (max_val <= INT32_MAX, 
        GrB_INVALID_VALUE, "Matrix values too large. Exiting to prevent overflow") ;
    GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n, n));
    GrB_Matrix_get_INT32(A, (int *) (&tcode), GrB_EL_TYPE_CODE);
    LG_ASSERT(tcode < 12 && tcode > 0, GrB_DOMAIN_MISMATCH);
    GRB_TRY (GrB_Matrix_assign (
            S, NULL, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL));
    if (sanitize)
    {
        // S = A+A'
        GRB_TRY (GrB_eWiseAdd (S, 0, 0, GrB_PLUS_UINT64, S, S, GrB_DESC_T1));
    }

    GRB_TRY (GrB_Matrix_new (&T, GrB_UINT64, n, n));
    GRB_TRY (GrB_Vector_new (&t, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&f, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&ramp, GrB_INT64, n + 1));
    GRB_TRY (GrB_Vector_new (&edge, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&cedge, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n));
    GRB_TRY (GrB_Vector_new (&index, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&parent_v, GrB_UINT64, n));

    // temporary arrays
    LG_TRY (LAGraph_Malloc ((void **) &V, n, sizeof (uint64_t), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SI, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SJ, 2*n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &SX, 2*n, sizeof (GrB_Index), msg)) ;

    // context arrays
    LG_TRY (LAGraph_Malloc ((void **) &context.data, 3 * n, sizeof (uint64_t), msg)) ;
    context.parent  = context.data;
    context.partner = context.data + n;
    context.weight  = context.data + 2 * n;

    // prepare vectors
    for (uint64_t i = 0; i < n; i++)
        context.parent[i] = i;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        f, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        f, NULL, NULL, GrB_ROWINDEX_INT64, f, (int64_t) 0, NULL));
    GRB_TRY (GrB_Vector_dup (&I, f));
    GRB_TRY (GrB_Vector_assign_UINT64 (
        ramp, NULL, NULL, (uint64_t) 0, GrB_ALL, n + 1, NULL));
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (
        ramp, NULL, NULL, GrB_ROWINDEX_INT64, ramp, (int64_t) 0, NULL));
    // semiring & monoid
    GrB_Index inf = ((uint64_t) INT_MAX << 32) ^ INT_MAX;
    GRB_TRY (GrB_BinaryOp_new (
        &comb, (GxB_binary_function) combine, 
        GrB_UINT64, GrB_UINT64, GrB_UINT64
    ));
    GRB_TRY (GrB_Semiring_new (&minComb, GrB_MIN_MONOID_UINT64, comb));
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
    GrB_Index nvals, ntuples = 0, num;
    bool diff = false;
    GRB_TRY (GrB_Matrix_nvals (&nvals, S));
    for (int iters = 1; nvals > 0; iters++)
    {
        // every vertex points to a root vertex at the beginning
        // edge[u] = u's minimum edge (weight and index are encoded together)
        GRB_TRY (GrB_assign (edge, 0, 0, inf, GrB_ALL, 0, 0));
        GRB_TRY (GrB_mxv (edge, 0, GrB_MIN_UINT64, minComb, S, f, 0));
        // cedge[u] = children's minimum edge  | if u is a root
        //          = (INT_MAX, u)             | otherwise
        GRB_TRY (GrB_assign (t, 0, 0, (uint64_t) INT_MAX, GrB_ALL, 0, 0));
        GRB_TRY (GrB_eWiseMult (cedge, 0, 0, comb, t, I, 0));
        LG_TRY (Reduce_assign (cedge, edge, context.parent, n, msg));
        // if (f[u] == u) f[u] := snd(cedge[u])  -- the index part of the edge
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, f, I, 0));
        GRB_TRY (GrB_apply (f, mask, GrB_SECOND_UINT64, snd, cedge, 0));
        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)
        GRB_TRY (GrB_Vector_extractTuples_UINT64 (NULL, V, &n, f));
        GRB_TRY (GrB_Vector_extract (t, NULL, NULL, f, V, n, NULL));
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, I, t, 0));
        GRB_TRY (GrB_assign (f, mask, GrB_MIN_UINT64, I, GrB_ALL, 0, 0));

        // five steps to generate the solution
        // 1. new roots (f[i] == i) revise their entries in cedge
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_EQ_UINT64, I, f, 0));
        GRB_TRY (GrB_assign (cedge, mask, 0, inf, GrB_ALL, 0, 0));

        // 2. every vertex tries to know whether one of its edges is selected
        GRB_TRY (GrB_extract (t, 0, 0, cedge, context.parent, n, 0));
        GRB_TRY (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, edge, t, 0));

        // 3. each root picks a vertex from its children to generate the solution
        GRB_TRY (GrB_assign (index, 0, 0, n, GrB_ALL, 0, 0));
        GRB_TRY (GrB_assign (index, mask, 0, I, GrB_ALL, 0, 0));
        GRB_TRY (GrB_assign (t, 0, 0, n, GrB_ALL, 0, 0));
        LG_TRY (Reduce_assign (t, index, context.parent, n, msg));
        GRB_TRY (GrB_extract (index, 0, 0, t, context.parent, n, 0));
        GRB_TRY (GrB_eWiseMult (mask ,0, 0, GrB_EQ_UINT64, I, index, 0));

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
        GRB_TRY (GrB_mxv (edge, mask, GrB_MIN_UINT64, minComb, T, I, 0));
        GRB_TRY (GrB_Vector_nvals (&num, edge));
        GRB_TRY (GrB_apply (t, 0, 0, snd, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (SI + ntuples, SJ + ntuples, &num, t));
        GRB_TRY (GrB_apply (t, 0, 0, fst, edge, 0));
        GRB_TRY (GrB_Vector_extractTuples (SI + ntuples, SX + ntuples, &num, t));
        GRB_TRY (GrB_Vector_clear (t));
        ntuples += num;

        // path halving until every vertex points on a root
        do {
            GRB_TRY (GrB_Vector_extractTuples_UINT64 (NULL, V, &n, f));
            GRB_TRY (GrB_Vector_extract (t, NULL, NULL, f, V, n, NULL));
            GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, f, t, 0));
            GrB_Vector temp = f;
            f = t;
            t = temp;
            temp = NULL;
            GRB_TRY (GrB_Vector_reduce_BOOL (&diff, NULL, GrB_LOR_MONOID_BOOL, mask, 0));
        } while (diff);

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
