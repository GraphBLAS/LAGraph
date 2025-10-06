//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_RPQMatrix.c: test cases for RPQ-matrix
// reachability algorithm
//------------------------------------------------------------------------------
//
// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Contributed by Rodion Suvorov, Semyon Grigoriev, St. Petersburg State
// University.

//------------------------------------------------------------------------------

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <acutest.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#define LEN 512
char msg[LAGRAPH_MSG_LEN] ;

static void load_matrix(GrB_Matrix *M, const char *name)
{
    char filename[LEN + 1] ;
    snprintf(filename, LEN, LG_DATA_DIR "%s", name) ;
    FILE *f = fopen(filename, "r") ;
    TEST_CHECK(f != NULL) ;
    OK(LAGraph_MMRead(M, f, msg)) ;
    fclose(f) ;
}

static char *matrix_to_str(GrB_Matrix M)
{
    GrB_Index nnz = 0 ;
    OK(GrB_Matrix_nvals(&nnz, M)) ;

    GrB_Index *I = NULL ;
    GrB_Index *J = NULL ;
    bool *X = NULL ;

    OK(LAGraph_Malloc((void **)&I, nnz, sizeof(GrB_Index), msg)) ;
    OK(LAGraph_Malloc((void **)&J, nnz, sizeof(GrB_Index), msg)) ;
    OK(LAGraph_Malloc((void **)&X, nnz, sizeof(bool), msg)) ;

    OK(GrB_Matrix_extractTuples_BOOL(I, J, X, &nnz, M)) ;

    size_t bufsize = 11 * nnz + 1 ;
    char *buf = NULL ;
    OK(LAGraph_Malloc((void **)&buf, bufsize, sizeof(char), msg)) ;
    buf[0] = '\0' ;

    for (size_t k = 0; k < nnz; k++)
    {
        char tmp[64] ;
        snprintf(tmp, sizeof(tmp), k == 0 ? "(%" PRIu64 ", %" PRIu64 ")" : " (%" PRIu64 ", %" PRIu64 ")", I[k] + 1, J[k] + 1) ;
        strcat(buf, tmp) ;
    }

    LAGraph_Free((void **)&I, msg) ;
    LAGraph_Free((void **)&J, msg) ;
    LAGraph_Free((void **)&X, msg) ;

    return buf ;
}

static void check_result_matrix(GrB_Matrix result, const char *expected)
{
    char *actual = matrix_to_str(result) ;
    TEST_CHECK(strcmp(actual, expected) == 0) ;
    if (strcmp(actual, expected) != 0)
    {
        TEST_MSG("Expected: %s\nActual:   %s", expected, actual) ;
    }
    LAGraph_Free((void **)&actual, msg) ;
}

//========================
// Tests with valid result
//========================

// a/b
void test_RPQMatrix_CONCAT(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A, B ;
    load_matrix(&A, "rpq_data/a.mtx") ;
    load_matrix(&B, "rpq_data/b.mtx") ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan concat = {
        .op = RPQ_MATRIX_OP_CONCAT,
        .lhs = &planA,
        .rhs = &planB} ;
    GrB_Index result_nnz ;
    OK(LAGraph_RPQMatrix(&result_nnz, &concat, msg)) ;

    const char *expected = "(1, 5) (1, 7) (2, 4) (2, 6) (3, 3) (7, 3)" ;
    check_result_matrix(concat.res_mat, expected) ;
    TEST_CHECK(result_nnz == 6) ;

    LAGraph_Finalize(msg) ;
}

// a|b
void test_RPQMatrix_LOR(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A, B ;
    load_matrix(&A, "rpq_data/a.mtx") ;
    load_matrix(&B, "rpq_data/b.mtx") ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan lor = {
        .op = RPQ_MATRIX_OP_LOR,
        .lhs = &planA,
        .rhs = &planB} ;
    GrB_Index result_nnz ;
    OK(LAGraph_RPQMatrix(&result_nnz, &lor, msg)) ;

    const char *expected = "(1, 2) (1, 3) (1, 7) (2, 4) (2, 5) (2, 7) (3, 6) (4, 4) (4, 6) "
                           "(5, 1) (5, 8) (6, 3) (7, 6)" ;
    check_result_matrix(lor.res_mat, expected) ;
    TEST_CHECK(result_nnz == 13) ;

    LAGraph_Finalize(msg) ;
}

// a*
void test_RPQMatrix_KLEENE(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A ;
    load_matrix(&A, "rpq_data/a.mtx") ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan kleene = {.op = RPQ_MATRIX_OP_KLEENE, .rhs = &planA} ;
    GrB_Index result_nnz ;
    OK(LAGraph_RPQMatrix(&result_nnz, &kleene, msg)) ;

    const char *expected =
        "(1, 1) (1, 2) (1, 4) (1, 6) (1, 7) (2, 2) (2, 4) (3, 3) "
        "(3, 6) (4, 4) (5, 5) (5, 8) (6, 6) (7, 6) (7, 7) (8, 8)" ;
    check_result_matrix(kleene.res_mat, expected) ;
    TEST_CHECK(result_nnz == 16) ;

    LAGraph_Finalize(msg) ;
}

// (a)*/b
void test_RPQMatrix_KLEENE_L(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A, B ;
    load_matrix(&A, "rpq_data/a.mtx") ;
    load_matrix(&B, "rpq_data/b.mtx") ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan kleene = {
        .op = RPQ_MATRIX_OP_KLEENE,
        .rhs = &planA} ;
    RPQMatrixPlan concat = {
        .op = RPQ_MATRIX_OP_CONCAT,
        .lhs = &kleene,
        .rhs = &planB} ;

    GrB_Index result_nnz ;
    OK(LAGraph_RPQMatrix(&result_nnz, &concat, msg)) ;
    const char *expected = "(1, 3) (1, 4) (1, 5) (1, 6) (1, 7) "
                           "(2, 4) (2, 5) (2, 6) (2, 7) (3, 3) "
                           "(4, 4) (4, 6) (5, 1) (6, 3) (7, 3)" ;
    check_result_matrix(concat.res_mat, expected) ;
    TEST_CHECK(result_nnz == 15) ;
    GrB_Info expected_nnz = result_nnz ;
    RPQMatrixPlan planA2 = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB2 = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan kleeneL = {
        .op = RPQ_MATRIX_OP_KLEENE_L,
        .lhs = &planA2,
        .rhs = &planB2} ;

    OK(LAGraph_RPQMatrix(&result_nnz, &kleeneL, msg)) ;

    expected = matrix_to_str(concat.res_mat) ;

    check_result_matrix(kleeneL.res_mat, expected) ;
    TEST_CHECK(result_nnz == expected_nnz) ;

    LAGraph_Finalize(msg) ;
}

// b/(a)*
void test_RPQMatrix_KLEENE_R(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A, B ;
    load_matrix(&A, "rpq_data/a.mtx") ;
    load_matrix(&B, "rpq_data/b.mtx") ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan kleene = {
        .op = RPQ_MATRIX_OP_KLEENE,
        .rhs = &planB} ;
    RPQMatrixPlan concat = {
        .op = RPQ_MATRIX_OP_CONCAT,
        .lhs = &planA,
        .rhs = &kleene} ;

    GrB_Index result_nnz;
    OK(LAGraph_RPQMatrix(&result_nnz, &concat, msg)) ;
    const char *expected = "(1, 1) (1, 2) (1, 3) (1, 5) (1, 7) "
                           "(2, 3) (2, 4) (2, 6) (3, 3) (3, 6) "
                           "(5, 8) (7, 3) (7, 6)" ;
    check_result_matrix(concat.res_mat, expected) ;
    TEST_CHECK(result_nnz == 13) ;
    GrB_Info expected_nnz = result_nnz ;
    RPQMatrixPlan planA2 = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB2 = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan kleeneR = {
        .op = RPQ_MATRIX_OP_KLEENE_R,
        .lhs = &planA2,
        .rhs = &planB2} ;

    OK(LAGraph_RPQMatrix(&result_nnz, &kleeneR, msg)) ;

    expected = matrix_to_str(concat.res_mat) ;
    check_result_matrix(kleeneR.res_mat, expected) ;
    TEST_CHECK(result_nnz == expected_nnz) ;

    LAGraph_Finalize(msg) ;
}

// c/(a|b)*
void test_RPQMatrix_Complex(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A, B, C ;
    load_matrix(&A, "rpq_data/a.mtx") ;
    load_matrix(&B, "rpq_data/b.mtx") ;
    load_matrix(&C, "rpq_data/c.mtx") ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan planC = {.op = RPQ_MATRIX_OP_LABEL, .mat = C} ;
    RPQMatrixPlan lor = {.op = RPQ_MATRIX_OP_LOR, .lhs = &planA, .rhs = &planB} ;
    RPQMatrixPlan kleene = {
        .op = RPQ_MATRIX_OP_KLEENE,
        .rhs = &lor} ;

    RPQMatrixPlan concat = {.op = RPQ_MATRIX_OP_CONCAT, .lhs = &planC, .rhs = &kleene} ;
    GrB_Index result_nnz ;
    OK(LAGraph_RPQMatrix(&result_nnz, &concat, msg)) ;

    const char *expected = "(1, 1) (1, 2) (1, 3) (1, 4) (1, 5) (1, 6) (1, 7) (1, 8) "
                           "(2, 3) (2, 6) (2, 8) "
                           "(3, 1) (3, 2) (3, 3) (3, 4) (3, 5) (3, 6) (3, 7) (3, 8) "
                           "(4, 1) (4, 2) (4, 3) (4, 4) (4, 5) (4, 6) (4, 7) (4, 8) "
                           "(5, 1) (5, 2) (5, 3) (5, 4) (5, 5) (5, 6) (5, 7) (5, 8) "
                           "(6, 3) (6, 4) (6, 6) "
                           "(7, 1) (7, 2) (7, 3) (7, 4) (7, 5) (7, 6) (7, 7) (7, 8) "
                           "(8, 3) (8, 6)" ;
    check_result_matrix(concat.res_mat, expected) ;
    TEST_CHECK(result_nnz == 48) ;
    LAGraph_Finalize(msg) ;
}

//==========================
// Tests with invalid result
//==========================

void test_RPQMatrix_Incorrect_size(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix M ;
    OK(GrB_Matrix_new(&M, GrB_BOOL, 3, 4)) ;

    RPQMatrixPlan plan = {.op = RPQ_MATRIX_OP_LABEL, .mat = M} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &plan, msg) ;
    TEST_CHECK(info == GrB_INVALID_VALUE) ;
    TEST_MSG("Expected GrB_INVALID_VALUE for non-square matrix") ;

    GrB_Matrix_free(&M) ;
    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Unequal_size(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A, B ;
    OK(GrB_Matrix_new(&A, GrB_BOOL, 3, 3)) ;
    OK(GrB_Matrix_new(&B, GrB_BOOL, 4, 4)) ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;
    RPQMatrixPlan planB = {.op = RPQ_MATRIX_OP_LABEL, .mat = B} ;
    RPQMatrixPlan concat = {.op = RPQ_MATRIX_OP_CONCAT, .lhs = &planA, .rhs = &planB} ;

    GrB_Index result_nnz ;
    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &concat, msg) ;

    TEST_CHECK(info == GrB_INVALID_VALUE) ;
    TEST_MSG("Expected GrB_INVALID_VALUE for unequal matrix sizes") ;

    GrB_Matrix_free(&A) ;
    GrB_Matrix_free(&B) ;
    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Null_Child_concat(void)
{
    LAGraph_Init(msg) ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL} ;
    RPQMatrixPlan concat = {.op = RPQ_MATRIX_OP_CONCAT, .lhs = &planA, .rhs = NULL} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &concat, msg) ;
    TEST_CHECK(info == GrB_NULL_POINTER || info == GrB_INVALID_VALUE) ;
    TEST_MSG("retval = %d (%s)", info, msg) ;

    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Null_Child_lor(void)
{
    LAGraph_Init(msg) ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL} ;
    RPQMatrixPlan lor = {.op = RPQ_MATRIX_OP_LOR, .lhs = &planA, .rhs = NULL} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &lor, msg) ;
    TEST_CHECK(info == GrB_NULL_POINTER || info == GrB_INVALID_VALUE) ;
    TEST_MSG("retval = %d (%s)", info, msg) ;

    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Null_Child_L_Kleene(void)
{
    LAGraph_Init(msg) ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL} ;
    RPQMatrixPlan kleeneL = {.op = RPQ_MATRIX_OP_KLEENE_L, .lhs = &planA, .rhs = NULL} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &kleeneL, msg) ;
    TEST_CHECK(info == GrB_NULL_POINTER || info == GrB_INVALID_VALUE) ;
    TEST_MSG("retval = %d (%s)", info, msg) ;

    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Null_Child_R_Kleene(void)
{
    LAGraph_Init(msg) ;

    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL} ;
    RPQMatrixPlan kleeneR = {.op = RPQ_MATRIX_OP_KLEENE_R, .lhs = &planA, .rhs = NULL} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &kleeneR, msg) ;
    TEST_CHECK(info == GrB_NULL_POINTER || info == GrB_INVALID_VALUE) ;
    TEST_MSG("retval = %d (%s)", info, msg) ;

    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Non_Null_Children_Kleene(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A ;
    OK(GrB_Matrix_new(&A, GrB_BOOL, 3, 3)) ;
    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;

    RPQMatrixPlan kleene = {.op = RPQ_MATRIX_OP_KLEENE, .lhs = &planA, .rhs = &planA} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &kleene, msg) ;
    TEST_CHECK(info == GrB_INVALID_VALUE) ;
    TEST_MSG("retval = %d (%s)", info, msg) ;

    GrB_Matrix_free(&A) ;
    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Left_Child_Kleene(void)
{
    LAGraph_Init(msg) ;

    GrB_Matrix A ;
    OK(GrB_Matrix_new(&A, GrB_BOOL, 3, 3)) ;
    RPQMatrixPlan planA = {.op = RPQ_MATRIX_OP_LABEL, .mat = A} ;

    RPQMatrixPlan kleene = {.op = RPQ_MATRIX_OP_KLEENE, .lhs = &planA, .rhs = NULL} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &kleene, msg) ;
    TEST_CHECK(info == GrB_NULL_POINTER) ;
    TEST_MSG("Expected GrB_NULL_POINTER for KLEENE with only left child") ;

    GrB_Matrix_free(&A) ;
    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Invalid_Op_type(void)
{
    LAGraph_Init(msg) ;

    RPQMatrixPlan plan = {.op = (RPQMatrixOp)52} ;
    GrB_Index result_nnz ;

    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, &plan, msg) ;
    TEST_CHECK(info == GrB_INVALID_VALUE) ;

    LAGraph_Finalize(msg) ;
}

void test_RPQMatrix_Null_root(void)
{
    LAGraph_Init(msg) ;

    GrB_Index result_nnz ;
    GrB_Info info = LAGraph_RPQMatrix(&result_nnz, NULL, msg) ;

    TEST_CHECK(info == GrB_NULL_POINTER) ;
    TEST_MSG("Expected GrB_NULL_POINTER for NULL plan root") ;

    LAGraph_Finalize(msg) ;
}


TEST_LIST = {
    {"RPQMatrix_CONCAT", test_RPQMatrix_CONCAT},
    {"RPQMatrix_LOR", test_RPQMatrix_LOR},
    {"RPQMatrix_KLEENE", test_RPQMatrix_KLEENE},
    {"RPQMatrix_KLEENE_L", test_RPQMatrix_KLEENE_L},
    {"RPQMatrix_KLEENE_R", test_RPQMatrix_KLEENE_R},
    {"RPQMatrix_Complex", test_RPQMatrix_Complex},
    {"RPQMatrix_Incorrect_size", test_RPQMatrix_Incorrect_size},
    {"RPQMatrix_Unequal_size", test_RPQMatrix_Unequal_size},
    {"RPQMatrix_Null_Child_concat", test_RPQMatrix_Null_Child_concat},
    {"RPQMatrix_Null_Child_lor", test_RPQMatrix_Null_Child_lor},
    {"RPQMatrix_Null_Child_L_Kleene", test_RPQMatrix_Null_Child_L_Kleene},
    {"RPQMatrix_Null_Child_R_Kleene", test_RPQMatrix_Null_Child_R_Kleene},
    {"RPQMatrix_Non_Null_Children_Kleene", test_RPQMatrix_Non_Null_Children_Kleene},
    {"RPQMatrix_Left_Child_Kleene", test_RPQMatrix_Left_Child_Kleene},
    {"RPQMatrix_Invalid_Op_type", test_RPQMatrix_Invalid_Op_type},
    {"RPQMatrix_Null_root", test_RPQMatrix_Null_root},
    {NULL, NULL}} ;
