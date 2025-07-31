#include <stdio.h>
#include <acutest.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg[LAGRAPH_MSG_LEN];
#define LEN 512

//****************************************************************************

void test_RPQMatrixConc(void)
{
    LAGraph_Init(msg);
    LAGraph_RpqMatrix_initialize();
    const char *nameA = "rpq_data/a.mtx";
    const char *nameB = "rpq_data/b.mtx";
    char filenameA [LEN+1] ;
    char filenameB [LEN+1] ;
    snprintf (filenameA, LEN, LG_DATA_DIR "%s", nameA) ;
    snprintf (filenameB, LEN, LG_DATA_DIR "%s", nameB) ;
    FILE *fA = fopen(filenameA, "r");
    FILE *fB = fopen(filenameB, "r");
    GrB_Matrix A, B;
    OK(LAGraph_MMRead(&A, fA, msg));
    OK(LAGraph_MMRead(&B, fB, msg));
    OK(fclose(fA));
    OK(fclose(fB));

    GrB_Index nvalsA, nvalsB;
    GrB_Matrix_nvals(&nvalsA,A);
    GrB_Matrix_nvals(&nvalsB,B);
    fprintf(stderr,"\nDEBUG: A:%lu and B:%lu\n",nvalsA,nvalsB);

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
    RpqMatrixPlan graphKleene = {
        .op = RPQ_MATRIX_OP_KLEENE,
        .lhs = NULL,
        .rhs = &graphConcat,
        .mat = NULL,
        .res_mat = NULL
    };    
    GrB_Index expected_nvasl = 14;
    GrB_Info res = LAGraph_RPQMatrix(&graphKleene,msg);
    GrB_Matrix result_matrix = graphKleene.res_mat;
    GrB_Index result;
    GrB_Matrix_nvals(&result,result_matrix);
    fprintf(stderr,"\nDEBUG: result: %lu",result);
    TEST_CHECK(result == expected_nvasl);
    LAGraph_Finalize(msg);
}

void test_RPQMatrixLor(void)
{
   LAGraph_Init(msg);
    LAGraph_RpqMatrix_initialize();
    const char *nameA = "rpq_data/a.mtx";
    const char *nameB = "rpq_data/b.mtx";
    char filenameA [LEN+1] ;
    char filenameB [LEN+1] ;
    snprintf (filenameA, LEN, LG_DATA_DIR "%s", nameA) ;
    snprintf (filenameB, LEN, LG_DATA_DIR "%s", nameB) ;
    FILE *fA = fopen(filenameA, "r");
    FILE *fB = fopen(filenameB, "r");
    GrB_Matrix A, B;
    OK(LAGraph_MMRead(&A, fA, msg));
    OK(LAGraph_MMRead(&B, fB, msg));
    OK(fclose(fA));
    OK(fclose(fB));

    GrB_Index nvalsA, nvalsB;
    GrB_Matrix_nvals(&nvalsA,A);
    GrB_Matrix_nvals(&nvalsB,B);
    fprintf(stderr,"\nDEBUG: A:%lu and B:%lu\n",nvalsA,nvalsB);
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
   RpqMatrixPlan graphKleene = {
        .op = RPQ_MATRIX_OP_KLEENE,
        .lhs = NULL,
        .rhs = &graphLor,
        .mat = NULL,
        .res_mat = NULL
    };    
    GrB_Index expected_nvasl = 35;
    GrB_Info res = LAGraph_RPQMatrix(&graphKleene,msg);
    GrB_Matrix result_matrix = graphKleene.res_mat;
    GrB_Index result;
    GrB_Matrix_nvals(&result,result_matrix);
    fprintf(stderr,"\nDEBUG: result: %lu",result);
    TEST_CHECK(result == expected_nvasl);
    LAGraph_Finalize(msg);
}

TEST_LIST = {
    // {"RPQMatrixKleene", test_RPQMatrixKleene},
    {"RPQMatrixConc", test_RPQMatrixConc},
    {"RPQMatrixLor", test_RPQMatrixLor},
    {NULL, NULL}};