//------------------------------------------------------------------------------
// LAGraph_scc.c
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

// TODO: not ready for src; uses global variables

/**
 * Code is based on the Min-Label algorithm described in the following paper:
 * D. Yan, J. Cheng, K. Xin, Y. Lu, W. Ng, Y. Bu, "Pregel Algorithms for Graph
 * Connectivity Problems with Performance Guarantees"
 * Proc. VLDB Endow. 7, 14 (October 2014), 1821–1832.
 * DOI: https://doi.org/10.14778/2733085.2733089
 **/

#define LG_FREE_ALL ;

#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

#if LAGRAPH_SUITESPARSE

//****************************************************************************
//arrays used in SelectOp
typedef struct 
{
    uint64_t *F, *B;
    bool *M;
} LG_SCC_Context;
#define SCCCONTEXT \
"typedef struct \n"             \
"{\n"                           \
"    uint64_t *F, *B;\n"    \
"    bool *M;\n"                \
"} LG_SCC_Context;\n"               


// LG_SCC_edge_removal:
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

void LG_SCC_edge_removal (bool *z, const void *x, GrB_Index i, GrB_Index j, const LG_SCC_Context *thunk) ;
void LG_SCC_edge_removal (bool *z, const void *x, GrB_Index i, GrB_Index j, const LG_SCC_Context *thunk)
{
    (*z) = (!thunk->M[i] && !thunk->M[j] 
        && thunk->F[i] == thunk->F[j] 
        && thunk->B[i] == thunk->B[j]) ;
}
#define EDGE_REMOVAL \
"void LG_SCC_edge_removal \n"                                                          \
"(bool *z, const void *x, GrB_Index i, GrB_Index j, const LG_SCC_Context *thunk)\n" \
"{\n"                                                                           \
"    (*z) = (!thunk->M[i] && !thunk->M[j] \n"                                   \
"        && thunk->F[i] == thunk->F[j] \n"                                      \
"        && thunk->B[i] == thunk->B[j]) ;\n"                                    \
"}\n"                                                                           

//****************************************************************************
// LG_SCC_trim_one: remove the edges connected to trivial SCCs
//  - A vertex is a trivial SCC if it has no incoming or outgoing edges.
//  - M[i] = i   | if vertex i is a trivial SCC
//    M[i] = n   | otherwise

void LG_SCC_trim_one (bool *z, const void *x, GrB_Index i, GrB_Index j, const LG_SCC_Context *thunk) ;
void LG_SCC_trim_one (bool *z, const void *x, GrB_Index i, GrB_Index j, const LG_SCC_Context *thunk)
{
    (*z) = (thunk->F[i] == thunk->F[j]) ;
}
#define TRIM_ONE \
"void LG_SCC_trim_one\n"                                                               \
"(bool *z, const void *x, GrB_Index i, GrB_Index j, const LG_SCC_Context *thunk)\n" \
"{\n"                                                                           \
"    (*z) = (thunk->F[i] == thunk->F[j]) ;\n"                                   \
"}\n"

//****************************************************************************
// label propagation
//  - label  : (input/output) labels
//  - mask   : (input) mask
//  - A      : (input) original matrix
//  - AT     : (input) transposed matrix
//  - n      : (input) number of vertices

#undef  LG_FREE_ALL
#define LG_FREE_ALL    \
{                      \
    GrB_free (&s) ;    \
    GrB_free (&t) ;    \
}

static GrB_Info propagate (GrB_Vector label, GrB_Vector mask,
        const GrB_Matrix A, const GrB_Matrix AT, GrB_Index n, char *msg)
{
    GrB_Info info;
    // semirings
    GrB_Vector s = NULL, t = NULL;
    GRB_TRY (GrB_Vector_new (&s, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&t, GrB_UINT64, n));
    GRB_TRY (GrB_assign (s, mask, 0, label, GrB_ALL, 0, 0));
    // GxB_fprint(s, GxB_SHORT, stdout);
    GRB_TRY (GrB_assign (t, 0, 0, label, GrB_ALL, 0, 0));
    GRB_TRY (GrB_wait(A, GrB_MATERIALIZE));

    bool active;
    while (true)
    {
        // GRB_TRY (GrB_mxv 
        //     (t, 0, GrB_MIN_UINT64, GrB_MIN_SECOND_SEMIRING_UINT64, AT, s, 0));
        GRB_TRY (GrB_vxm (t, 0, GrB_MIN_UINT64,
                                 GrB_MIN_FIRST_SEMIRING_UINT64, s, A, 0));
        GRB_TRY (GrB_eWiseMult (mask, 0, 0, GrB_NE_UINT64, t, label, 0));
        GRB_TRY (GrB_assign (label, NULL, NULL, t, GrB_ALL, n, NULL));
        GRB_TRY (GrB_reduce (&active, 0, GrB_LOR_MONOID_BOOL, mask, 0));
        if (!active) break;
        GRB_TRY (GrB_Vector_clear(s));
        GRB_TRY (GrB_assign (s, mask, 0, label, GrB_ALL, 0, 0));
        GRB_TRY (GrB_wait(s, GrB_MATERIALIZE));
    }

    LG_FREE_ALL ;
    return GrB_SUCCESS;
}
//****************************************************************************

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
    LAGraph_Free ((void **) &contx.F, msg); \
    LAGraph_Free ((void **) &contx.B, msg); \
    LAGraph_Free ((void **) &contx.M, msg); \
    GrB_free (&ind);                        \
    GrB_free (&inf);                        \
    GrB_free (&f);                          \
    GrB_free (&b);                          \
    GrB_free (&D);                          \
    GrB_free (&x);                          \
    GrB_free (&mask);                       \
    GrB_free (&m2);                         \
    GrB_free (&FW);                         \
    GrB_free (&BW);                         \
    GrB_free (&sel1);                       \
    GrB_free (&sel2);                       \
    GrB_free (&contx_type);                 \
    GrB_free (&scc);

#endif

//****************************************************************************
int LAGraph_scc
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    LG_CLEAR_MSG ;
    LG_SCC_Context contx = {NULL, NULL, NULL};
    GrB_Info info = GrB_SUCCESS;
    GrB_Type contx_type = NULL;
    GrB_Type type_F = NULL, type_B = NULL, type_M = NULL;
    int hand_F = GrB_DEFAULT, hand_B = GrB_DEFAULT, hand_M = GrB_DEFAULT;
    uint64_t n_F = 0, n_B = 0, n_M = 0, size_F = 0, size_B = 0, size_M = 0;
    GrB_Vector scc = NULL ;
    GrB_Vector ind = NULL ;
    GrB_Vector inf = NULL ;
    GrB_Vector x = NULL ;
    GrB_Vector f = NULL, b = NULL, mask = NULL, m2 = NULL;
    GrB_IndexUnaryOp sel1 = NULL, sel2 = NULL ;
    GrB_Monoid Add = NULL ;
    GrB_Matrix FW = NULL, BW = NULL, D = NULL;
    LG_ASSERT(result != NULL, GrB_NULL_POINTER);
    LG_ASSERT(A != NULL, GrB_NULL_POINTER);

    GrB_Index n, ncols, nvals;
    GRB_TRY (GrB_Matrix_nrows (&n, A));
    GRB_TRY (GrB_Matrix_ncols (&ncols, A));
    LG_ASSERT(n == ncols, GrB_DIMENSION_MISMATCH);
    
    #if !LG_SUITESPARSE_GRAPHBLAS_V10
    LG_TRY (LAGraph_Malloc ((void **) &contx.F, n, sizeof (uint64_t), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &contx.B, n, sizeof (uint64_t), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &contx.M, n, sizeof (bool), msg)) ;
    #endif
    // scc: the SCC identifier for each vertex
    // scc[u] == n: not assigned yet
    GRB_TRY (GrB_Vector_new (&scc, GrB_UINT64, n));
    // vector of indices: ind[i] == i
    GRB_TRY (GrB_Vector_new (&ind, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_assign_UINT64 (
        ind, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (
        ind, NULL, NULL, GrB_ROWINDEX_INT64, ind, 0, NULL)) ;
    // vector of infinite value: inf[i] == n
    GRB_TRY (GrB_Vector_new (&inf, GrB_UINT64, n));
    GRB_TRY (GrB_assign (inf, NULL, NULL, n, GrB_ALL, 0, NULL));
    // other vectors
    GRB_TRY (GrB_Vector_new (&f, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&b, GrB_UINT64, n));
    GRB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n));
    GRB_TRY (GrB_Vector_new (&m2, GrB_BOOL, n));
    GRB_TRY (GrB_Vector_new (&x, GrB_BOOL, n));
    GRB_TRY (GxB_Type_new (
        &contx_type, sizeof(LG_SCC_Context), "LG_SCC_Context", SCCCONTEXT)) ;

    GRB_TRY (GxB_IndexUnaryOp_new (
        &sel1, (GxB_index_unary_function) LG_SCC_trim_one, 
        GrB_BOOL, GrB_UINT64, contx_type, 
        // NULL, NULL
        "LG_SCC_trim_one", TRIM_ONE
    ));
    GRB_TRY (GxB_IndexUnaryOp_new (
        &sel2, (GxB_index_unary_function) LG_SCC_edge_removal, 
        GrB_BOOL, GrB_UINT64, contx_type,
        // NULL, NULL
        "LG_SCC_edge_removal", EDGE_REMOVAL
    ));

    // store the graph in both directions (forward / backward)
    GRB_TRY (GrB_Matrix_new (&FW, GrB_BOOL, n, n));
    GRB_TRY (GrB_Matrix_new (&BW, GrB_BOOL, n, n));
    GRB_TRY (GrB_Matrix_assign_BOOL(
        FW, A, NULL, true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
    GRB_TRY (GrB_transpose (BW, NULL, NULL, FW, NULL));     // BW = FW'
       
    // check format
    int32_t A_format, AT_format;
    GRB_TRY (GrB_get (FW, &A_format , GrB_STORAGE_ORIENTATION_HINT));
    GRB_TRY (GrB_get (BW, &AT_format, GrB_STORAGE_ORIENTATION_HINT));

    bool is_csr = (A_format == GrB_ROWMAJOR && AT_format == GrB_ROWMAJOR);
    LG_ASSERT (is_csr, GrB_INVALID_VALUE) ;

    // remove trivial SCCs
    GRB_TRY (GrB_Vector_assign_BOOL (x, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_mxv (m2, NULL, NULL, GxB_ANY_PAIR_BOOL, FW, x, NULL)) ;
    GRB_TRY (GrB_mxv (mask, m2, NULL, GxB_ANY_PAIR_BOOL, BW, x, GrB_DESC_S)) ;
    GRB_TRY (GrB_Vector_nvals (&nvals, mask));

    GRB_TRY (GrB_assign (scc, NULL, NULL, ind, GrB_ALL, 0, NULL));
    GRB_TRY (GrB_assign (scc, mask, NULL, n, GrB_ALL, 0, NULL));

    if (nvals < n)
    {
        // No reason for context. 
        #if LG_SUITESPARSE_GRAPHBLAS_V10
        GRB_TRY(GxB_Vector_unload(
            scc, (void **) &contx.F, &type_F, &n_F, &size_F, &hand_F, NULL)) ;
        #else
        GRB_TRY (GrB_Vector_extractTuples_UINT64 (NULL, contx.F, &n, scc));
        #endif
        GRB_TRY (GrB_Matrix_select_UDT (
            FW, NULL, NULL, sel1, FW, &contx, NULL)) ;
        GRB_TRY (GrB_Matrix_select_UDT (
            BW, NULL, NULL, sel1, BW, &contx, NULL)) ;
        #if LG_SUITESPARSE_GRAPHBLAS_V10
        GRB_TRY(GxB_Vector_load(
            scc, (void **) &contx.F, type_F, n_F, size_F, hand_F, NULL)) ;
        #endif
    }

    GRB_TRY (GrB_Matrix_nvals (&nvals, FW));
    while (nvals > 0)
    {
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64 (
            mask, NULL, NULL, GrB_EQ_UINT64, scc, n, NULL));
        GRB_TRY (GrB_assign (f, NULL, NULL, ind, GrB_ALL, 0, NULL));
        LG_TRY (propagate (f, mask, FW, BW, n, msg));

        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, f, ind, NULL));
        GRB_TRY (GrB_Vector_assign_UINT64 (
            b, NULL, NULL, n, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_assign (b, mask, NULL, ind, GrB_ALL, 0, NULL));
        LG_TRY (propagate (b, mask, BW, FW, n, msg));

        GRB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, f, b, NULL));
        GRB_TRY (GrB_assign (scc, mask, GrB_MIN_UINT64, f, GrB_ALL, 0, NULL));

        #if LG_SUITESPARSE_GRAPHBLAS_V10
            GRB_TRY(GxB_Vector_unload(
                f, (void **) &contx.F, &type_F, &n_F, &size_F, &hand_F, NULL)) ;
            GRB_TRY(GxB_Vector_unload(
                b, (void **) &contx.B, &type_B, &n_B, &size_B, &hand_B, NULL)) ;
            GRB_TRY(GxB_Vector_unload(
                mask, (void **) &contx.M, &type_M, &n_M, &size_M, &hand_M, NULL
            )) ;
        #else
            GRB_TRY (GrB_Vector_extractTuples_UINT64 (NULL, contx.F, &n, f));
            GRB_TRY (GrB_Vector_extractTuples_UINT64 (NULL, contx.B, &n, b));
            GRB_TRY (GrB_Vector_extractTuples_BOOL (NULL, contx.M, &n, mask));
        #endif

        GRB_TRY (GrB_Matrix_select_UDT (
            FW, NULL, NULL, sel2, FW, &contx, NULL)) ;
        GRB_TRY (GrB_Matrix_select_UDT (
            BW, NULL, NULL, sel2, BW, &contx, NULL)) ;
        #if LG_SUITESPARSE_GRAPHBLAS_V10
            GRB_TRY(GxB_Vector_load(
                f, (void **) &contx.F, type_F, n_F, size_F, hand_F, NULL)) ;
            GRB_TRY(GxB_Vector_load(
                b, (void **) &contx.B, type_B, n_B, size_B, hand_B, NULL)) ;
            GRB_TRY(GxB_Vector_load(
                mask, (void **) &contx.M, type_M, n_M, size_M, hand_M, NULL
            )) ;
        #endif
        GRB_TRY (GrB_Matrix_nvals (&nvals, FW));
    }
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64 (
            mask, NULL, NULL, GrB_EQ_UINT64, scc, n, NULL));
    GRB_TRY (GrB_assign (scc, mask, NULL, ind, GrB_ALL, 0, NULL));

    *result = scc;
    scc = NULL;

    LG_FREE_ALL ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
