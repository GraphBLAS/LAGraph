#include <stdio.h>
#include <acutest.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg[LAGRAPH_MSG_LEN];

//****************************************************************************
void test_RPQMatrixKleene(void)
{
    LAGraph_Init(msg);
    const char *nameA = "rpq_data/a.mtx";
    FILE *fA = fopen(nameA, "r");
    TEST_CHECK(fA != NULL);
    GrB_Matrix A;
    OK(LAGraph_MMRead(&A, fA, msg));
    OK(fclose(fA));
    RpqMatrixPlan graphA = {
        .op = RPQ_MATRIX_OP_LABEL,
        .lhs = NULL,
        .rhs = NULL,
        .mat = A,
        .res_mat = NULL
    };
    RpqMatrixPlan graphKleene = {
        .op = RPQ_MATRIX_OP_KLEENE,
        .lhs = NULL,
        .rhs = &graphA,
        .mat = NULL,
        .res_mat = NULL
    };
    GrB_Info res = LAGraph_RpqMatrix(&graphKleene,msg);
    GrB_Matrix result_matrix = graphKleene.res_mat;
    GrB_Index result;
    GrB_Matrix_nvals(&result,result_matrix);
    LAGraph_Finalize(msg);
}

void test_RPQMatrixConc(void)
{
    LAGraph_Init(msg);
    const char *nameA = "rpq_data/a.mtx";
    const char *nameB = "rpq_data/b.mtx";
    FILE *fA = fopen(nameA, "r");
    FILE *fB = fopen(nameB, "r");
    TEST_CHECK(fA != NULL);
    TEST_CHECK(fB != NULL);
    GrB_Matrix A, B;
    OK(LAGraph_MMRead(&A, fA, msg));
    OK(LAGraph_MMRead(&B, fB, msg));
    OK(fclose(fA));
    OK(fclose(fB));
    RpqMatrixPlan graphA = {
        .op = RPQ_MATRIX_OP_LABEL,
        .lhs = NULL,
        .rhs = NULL,
        .mat = A,
        .res_mat = NULL
    };
    RpqMatrixPlan graphB = {
        .op = RPQ_MATRIX_OP_LABEL,
        .lhs = NULL,
        .rhs = NULL,
        .mat = B,
        .res_mat = NULL
    };
    RpqMatrixPlan graphConcat = {
        .op = RPQ_MATRIX_OP_CONCAT,
        .lhs = &graphA,
        .rhs = &graphB,
        .mat = NULL,
        .res_mat = NULL
    };
    GrB_Info res = LAGraph_RpqMatrix(&graphConcat,msg);
    GrB_Matrix result_matrix = graphConcat.res_mat;
    GrB_Index result;
    GrB_Matrix_nvals(&result,result_matrix);
    LAGraph_Finalize(msg);
}

void test_RPQMatrixLor(void)
{
    LAGraph_Init(msg);
    const char *nameA = "rpq_data/a.mtx";
    const char *nameB = "rpq_data/b.mtx";
    FILE *fA = fopen(nameA, "r");
    FILE *fB = fopen(nameB, "r");
    TEST_CHECK(fA != NULL);
    TEST_CHECK(fB != NULL);
    GrB_Matrix A, B;
    OK(LAGraph_MMRead(&A, fA, msg));
    OK(LAGraph_MMRead(&A, fB, msg));
    OK(fclose(fA));
    OK(fclose(fB));
    RpqMatrixPlan graphA = {
        .op = RPQ_MATRIX_OP_LABEL,
        .lhs = NULL,
        .rhs = NULL,
        .mat = A,
        .res_mat = NULL
    };
    RpqMatrixPlan graphB = {
        .op = RPQ_MATRIX_OP_LABEL,
        .lhs = NULL,
        .rhs = NULL,
        .mat = B,
        .res_mat = NULL
    };
    RpqMatrixPlan graphLor = {
        .op = RPQ_MATRIX_OP_LOR,
        .lhs = &graphA,
        .rhs = &graphB,
        .mat = NULL,
        .res_mat = NULL
    };
        GrB_Info res = LAGraph_RpqMatrix(&graphLor,msg);
    GrB_Matrix result_matrix = graphLor.res_mat;
    GrB_Index result;
    GrB_Matrix_nvals(&result,result_matrix);
    LAGraph_Finalize(msg);
}

TEST_LIST = {
    {"RPQMatrixKleene", test_RPQMatrixKleene},
    {"RPQMatrixConc", test_RPQMatrixConc},
    {"RPQMatrixLor", test_RPQMatrixLor},
    {NULL, NULL}};