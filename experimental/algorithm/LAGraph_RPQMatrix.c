//------------------------------------------------------------------------------
// LAGraph_RPQMatrix: regular path query algortithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// Contributed by Rodion Suvorov, Semyon Grigoriev, St. Petersburg State
// University.
//
//------------------------------------------------------------------------------

// Code is based on the algorithm described in the following paper:
//  * Diego Arroyuelo, Adrián Gómez-Brandón & Gonzalo Navarro "Evaluating
//    regular path queries on compressed adjacency matrices"
//  * URL: https://link.springer.com/article/10.1007/s00778-024-00885-6

//------------------------------------------------------------------------------
// LAGraph_RPQMatrix: regular path query algortithm
//
// For an edge-labelled directed graph the algorithm computes the nubmer of
// nonzero elements in its reachability matrix.
// The reachability matrix created by following rules:
// * A[i,j] = True if node with index j is reachable from node with index i
//   and concatenation of labels over path between these two labels is a word
//   from specified regular language.
// * A[i,j] = False in other cases.
//
// The algorithm is based on the idea of ​​considering a regular constraint as
// an abstract syntax tree, the leaves of which are matrices of adjacency matrix
// decomposition of the graph, and the internal nodes are the operations of
// conjunction, concatenation, etc.
//
// Example of adjacency matrix decomposition:
//
// Graph:
// (0) --[a]-> (1)
//  |           ^
// [b]    [c]--/
//  |  --/
//  v /
// (2) --[b]-> (3)
//
// Adjacency matrix decomposition of this graph consists of:
// * Adjacency matrix for the label a:
//       0   1   2   3
//   0 |   | T |   |   |
//   1 |   |   |   |   |
//   2 |   |   |   |   |
//   3 |   |   |   |   |
// * Adjacency matrix for the label b:
//       0   1   2   3
//   0 |   |   | T |   |
//   1 |   |   |   |   |
//   2 |   |   |   | T |
//   3 |   |   |   |   |
// * Adjacency matrix for the label c:
//       0   1   2   3
//   0 |   |   |   |   |
//   1 |   |   |   |   |
//   2 |   | T |   |   |
//   3 |   |   |   |   |
//
// The algorithm recursively starts from the root of the given tree and
// performs the operations corresponding to each node on the children of that
// node. As a result of the algorithm's execution, the reachability
// matrix will be stored at the root.
//
// Example of regular expression and its corresponding AST:
//
// Regular expression:
// a/(b|c)*
//
// Abstract syntax tree:
//    ┌─┐
//    │/| (3)
//    └┬┘
// ┌─┬─┴─┬─┐
// │a│   │*│ (2)
// └─┘   └┬┘
//       ┌┴┐
//       │|│ (1)
//       └┬┘
//    ┌─┬─┴─┬─┐
//    │b│   │c│
//    └─┘   └─┘
// The numbers next to the graph nodes show the order in which operations are
// executed. For the decomposition and AST specified above, the resulting
// matrix will have the following structure (Note, that * represents the
// reflexive-transitive closure):
//
//      0   1   2   3
//  0 |   | T |   |   |
//  1 |   |   |   |   |
//  2 |   |   |   |   |
//  3 |   |   |   |   |
//
// So for this example LAGraph_RPQMatrix will return 1.
//
// Full description available at:
//   https://arxiv.org/pdf/2307.14930

#define LG_FREE_WORK \
    {                \
    }

#define LG_FREE_ALL   \
    {                 \
        LG_FREE_WORK; \
    }

#include "LG_internal.h"
#include "LAGraphX.h"
#include <assert.h>

#define OK(s)                                               \
    {                                                       \
        GrB_Info info = s;                                  \
        if (!(info == GrB_SUCCESS))                         \
        {                                                   \
            printf("GraphBLAS error: %d\n", info);          \
            fprintf(stderr, "GraphBLAS error: %d\n", info); \
        }                                                   \
    }

#include <stdbool.h>
#include <stdio.h>

GrB_Info LAGraph_RPQMatrix_check(RPQMatrixPlan *plan, GrB_Index *dimension, char *msg)
{
    if (plan == NULL)
    {
        return GrB_SUCCESS;
    }
    if (plan->op == RPQ_MATRIX_OP_LABEL)
    {
        GrB_Index nrows, ncols;
        OK(GrB_Matrix_nrows(&nrows, plan->mat));
        OK(GrB_Matrix_ncols(&ncols, plan->mat));
        if (*dimension == -1)
        {
            LG_ASSERT_MSG(nrows == ncols, GrB_INVALID_VALUE,
                          "all the matrices in the graph adjacency matrix decomposition "
                          "should have the same dimensions and be square");
            *dimension = ncols;
        }
        else
        {
            LG_ASSERT_MSG(nrows == *dimension || ncols == *dimension, GrB_INVALID_VALUE,
                          "all the matrices in the graph adjacency matrix decomposition "
                          "should have the same dimensions and be square");
        }

        return GrB_SUCCESS;
    }
    GrB_Info lstatus = LAGraph_RPQMatrix_check(plan->lhs, dimension, msg);
    GrB_Info rstatus = LAGraph_RPQMatrix_check(plan->rhs, dimension, msg);
    if (rstatus || lstatus)
    {
        return GrB_INVALID_VALUE;
    }
    return GrB_SUCCESS;
}

static GrB_Semiring sr;
static GrB_Monoid op;

GrB_Info LAGraph_RPQMatrix_solver(RPQMatrixPlan *plan, char *msg);

static GrB_Info LAGraph_RPQMatrixLor(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_LOR, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RPQMatrixPlan *lhs = plan->lhs;
    RPQMatrixPlan *rhs = plan->rhs;

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix_solver(lhs, msg));
    OK(LAGraph_RPQMatrix_solver(rhs, msg));

    GrB_Matrix lhs_mat = lhs->res_mat;
    GrB_Matrix rhs_mat = rhs->res_mat;

    GrB_Index width, height;
    GrB_Matrix_nrows(&height, rhs_mat);
    GrB_Matrix_ncols(&width, lhs_mat);
    GrB_Matrix res;
    GrB_Matrix_new(&res, GrB_BOOL, height, width);
    GRB_TRY(GrB_eWiseAdd(res, GrB_NULL, GrB_NULL,
                         op, lhs_mat, rhs_mat, GrB_DESC_R));
    plan->res_mat = res;

    return (GrB_SUCCESS);
}

static GrB_Info LAGraph_RPQMatrixConcat(RPQMatrixPlan *plan, char *msg)
{

    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_CONCAT, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RPQMatrixPlan *lhs = plan->lhs;
    RPQMatrixPlan *rhs = plan->rhs;

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix_solver(lhs, msg));
    OK(LAGraph_RPQMatrix_solver(rhs, msg));

    GrB_Matrix lhs_mat = lhs->res_mat;
    GrB_Matrix rhs_mat = rhs->res_mat;

    GrB_Index width, height;
    GrB_Matrix_nrows(&height, rhs_mat);
    GrB_Matrix_ncols(&width, lhs_mat);
    GrB_Matrix res;
    GrB_Matrix_new(&res, GrB_BOOL, height, width);
    GRB_TRY(GrB_mxm(res, GrB_NULL, GrB_NULL,
                    sr, lhs_mat, rhs_mat, GrB_DESC_R));
    // GrB_mxm(plan->res_mat, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, lhs_mat, rhs_mat, GrB_NULL);
    plan->res_mat = res;

    return (GrB_SUCCESS);
}

static GrB_Info LAGraph_RPQMatrixKleene(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_KLEENE, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RPQMatrixPlan *lhs = plan->lhs;
    RPQMatrixPlan *rhs = plan->rhs;

    // Kleene star should have one child. Always right.
    LG_ASSERT(lhs == NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix_solver(rhs, msg));

    GrB_Matrix B = rhs->res_mat;

    // Creating identity matrix.
    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, B));
    GrB_Matrix I;
    GRB_TRY(GrB_Matrix_new(&I, GrB_BOOL, n, n));

    GrB_Vector v;
    GRB_TRY(GrB_Vector_new(&v, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Matrix_diag(&I, v, 0));

    GRB_TRY(GrB_Vector_free(&v));

    // B + I
    GrB_Matrix BPI;
    GRB_TRY(GrB_Matrix_new(&BPI, GrB_BOOL, n, n));
    GRB_TRY(GrB_eWiseAdd(BPI, GrB_NULL, GrB_NULL,
                         op, B, I, GrB_DESC_R));
    // S <- I
    GrB_Matrix S;
    GRB_TRY(GrB_Matrix_dup(&S, I));

    bool changed = true;
    GrB_Index nnz_S = n, nnz_Sold = 0;

    while (changed)
    {
        // S <- S x (B + I)
        GRB_TRY(GrB_mxm(S, GrB_NULL, GrB_NULL,
                        sr, S, BPI, GrB_DESC_R));

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S));
        if (nnz_S != nnz_Sold)
        {
            changed = true;
            nnz_Sold = nnz_S;
        }
        else
        {
            changed = false;
        }
    }
    plan->res_mat = S;

    GRB_TRY(GrB_Matrix_free(&I));
    GRB_TRY(GrB_Matrix_free(&BPI));
    return (GrB_SUCCESS);
}

// this function need to handle special case where some optimization
// are available.
//
// consider following AST:
//    ┌─┐
//    │/|
//    └┬┘
// ┌─┬─┴─┬─┐
// │*│   │b│
// └┬┘   └─┘
// ┌┴┐
// │a│
// └─┘
// If matrix B is sparse and A is dense, then instead of naive
// way:
//
// (I + A + A x A + ...) x B
//
// we can do:
//
// (B + A x B + A x A x B + ...)
//
// and AST should be rewritten in the following way:
//   ┌───┐
//   │L^*│
//   └─┬─┘
// ┌─┬─┴─┬─┐
// │a│   │b│
// └─┘   └─┘
static GrB_Info LAGraph_RPQMatrixKleene_L(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_KLEENE, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RPQMatrixPlan *lhs = plan->lhs; // A
    RPQMatrixPlan *rhs = plan->rhs; // B

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix_solver(lhs, msg));
    OK(LAGraph_RPQMatrix_solver(rhs, msg));

    GrB_Matrix A = lhs->res_mat;
    GrB_Matrix B = rhs->res_mat;

    // Creating identity matrix.
    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GrB_Matrix I;
    GRB_TRY(GrB_Matrix_new(&I, GrB_BOOL, n, n));

    GrB_Vector v;
    GRB_TRY(GrB_Vector_new(&v, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Matrix_diag(&I, v, 0));

    GRB_TRY(GrB_Vector_free(&v));

    // B + I
    GrB_Matrix API;
    GRB_TRY(GrB_Matrix_new(&API, GrB_BOOL, n, n));
    GRB_TRY(GrB_eWiseAdd(API, GrB_NULL, GrB_NULL,
                         op, A, I, GrB_DESC_R));

    // S <- B
    GrB_Matrix S;
    GRB_TRY(GrB_Matrix_dup(&S, B));

    bool changed = true;
    GrB_Index nnz_S = 0, nnz_Sold = 0;

    while (changed)
    {
        // T <- (A + I) x S
        GRB_TRY(GrB_mxm(S, NULL, NULL, sr, API, S, NULL));

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S));
        if (nnz_S != nnz_Sold)
        {
            changed = true;
            nnz_Sold = nnz_S;
        }
        else
        {
            changed = false;
        }
    }

    plan->res_mat = S;
    GRB_TRY(GrB_Matrix_free(&I));
    GRB_TRY(GrB_Matrix_free(&API));
    return GrB_SUCCESS;
}

// this function need to handle special case where some optimization
// are available.
// consider following AST:
//    ┌─┐
//    │/|
//    └┬┘
// ┌─┬─┴─┬─┐
// │a│   │*│
// └─┘   └┬┘
//       ┌┴┐
//       │b│
//       └─┘
// If matrix A is sparse and B is dense, then instead of naive
// way:
//
// A x (I + B + B x B + ...)
//
// we can do:
//
// (A + A x B + A x B x B + ...)
//
// and AST should be rewritten in the following way:
//   ┌───┐
//   │R^*│
//   └─┬─┘
// ┌─┬─┴─┬─┐
// │a│   │b│
// └─┘   └─┘

static GrB_Info LAGraph_RPQMatrixKleene_R(RPQMatrixPlan *plan, char *msg)
{
    LG_ASSERT(plan != NULL, GrB_NULL_POINTER);
    LG_ASSERT(plan->op == RPQ_MATRIX_OP_KLEENE, GrB_INVALID_VALUE);
    LG_ASSERT(plan->res_mat == NULL, GrB_INVALID_VALUE);

    RPQMatrixPlan *lhs = plan->lhs; // A
    RPQMatrixPlan *rhs = plan->rhs; // B

    LG_ASSERT(lhs != NULL, GrB_NULL_POINTER);
    LG_ASSERT(rhs != NULL, GrB_NULL_POINTER);

    OK(LAGraph_RPQMatrix_solver(lhs, msg));
    OK(LAGraph_RPQMatrix_solver(rhs, msg));

    GrB_Matrix A = lhs->res_mat;
    GrB_Matrix B = rhs->res_mat;

    // Creating identity matrix.
    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, B));
    GrB_Matrix I;
    GRB_TRY(GrB_Matrix_new(&I, GrB_BOOL, n, n));

    GrB_Vector v;
    GRB_TRY(GrB_Vector_new(&v, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Matrix_diag(&I, v, 0));

    GRB_TRY(GrB_Vector_free(&v));

    // B + I
    GrB_Matrix BPI;
    GRB_TRY(GrB_Matrix_new(&BPI, GrB_BOOL, n, n));
    GRB_TRY(GrB_eWiseAdd(BPI, GrB_NULL, GrB_NULL,
                         op, B, I, GrB_DESC_R));

    // S <- A
    GrB_Matrix S;
    GRB_TRY(GrB_Matrix_dup(&S, A));

    bool changed = true;
    GrB_Index nnz_S = 0, nnz_Sold = 0;

    while (changed)
    {
        // S <- S x (B + I)
        GRB_TRY(GrB_mxm(S, NULL, NULL, sr, S, BPI, NULL));

        GRB_TRY(GrB_Matrix_nvals(&nnz_S, S));
        if (nnz_S != nnz_Sold)
        {
            changed = true;
            nnz_Sold = nnz_S;
        }
        else
        {
            changed = false;
        }
    }

    plan->res_mat = S;
    GRB_TRY(GrB_Matrix_free(&I));
    GRB_TRY(GrB_Matrix_free(&BPI));
    return GrB_SUCCESS;
}

GrB_Info LAGraph_RPQMatrix_solver(RPQMatrixPlan *plan, char *msg)
{
    if (plan->res_mat != NULL)
    {
        return (GrB_SUCCESS);
    }

    switch (plan->op)
    {
    case RPQ_MATRIX_OP_LABEL:
        LG_ASSERT_MSG(plan->lhs == NULL && plan->rhs == NULL,
                      GrB_INVALID_VALUE, "label node should not have any children nodes");
        plan->res_mat = plan->mat;
        return (GrB_SUCCESS);
    case RPQ_MATRIX_OP_LOR:
        return LAGraph_RPQMatrixLor(plan, msg);
    case RPQ_MATRIX_OP_CONCAT:
        return LAGraph_RPQMatrixConcat(plan, msg);
    case RPQ_MATRIX_OP_KLEENE:
        return LAGraph_RPQMatrixKleene(plan, msg);
    case RPQ_MATRIX_OP_KLEENE_L:
        return LAGraph_RPQMatrixKleene_L(plan, msg);
    case RPQ_MATRIX_OP_KLEENE_R:
        return LAGraph_RPQMatrixKleene_R(plan, msg);
    default:
        LG_ASSERT_MSG(false, GrB_INVALID_VALUE, "invalid graph node type");
    }
    return (GrB_SUCCESS);
}

GrB_Info LAGraph_RPQMatrix_initialize()
{
    sr = GrB_LOR_LAND_SEMIRING_BOOL;
    op = GxB_LOR_BOOL_MONOID;
}

GrB_Info LAGrah_RPQMatrix(
    // output:
    GrB_Index *nnz, // number of nonzero values in
                    // result reachability matrix

    // input:
    RPQMatrixPlan *plan, // root of abstarct syntax tree of
                         // regular expression
    char *msg            // LAGraph output message
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;
    LG_ASSERT(plan == NULL, GrB_NULL_POINTER);
    GrB_Index dimension = -1;
    GrB_Info info = LAGraph_RPQMatrix_check(plan, &dimension, msg);
    LG_ASSERT_MSG(info == GrB_SUCCESS, info, msg);

    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------

    LAGraph_RPQMatrix_initialize();

    //--------------------------------------------------------------------------
    // run solver
    //--------------------------------------------------------------------------

    info = LAGraph_RPQMatrix_solver(plan, msg);
    LG_ASSERT_MSG(info == GrB_SUCCESS, info, msg);
    GrB_Matrinx_nvals(nnz, plan->res_mat);
    return GrB_SUCCESS;
}