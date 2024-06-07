#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

void test_MCM()
{
    LAGraph_Init(msg);
    GrB_Vector mateC = NULL;
    OK(GrB_Vector_new(&mateC, GrB_UINT64, 5));
    GrB_Index Ilist[2] = {3, 4};
    OK(GrB_Vector_assign_UINT64(mateC, NULL, NULL, 1, Ilist, 2, NULL));

    GrB_Index R[9] = {0, 0, 1, 2, 2, 3, 3, 4, 4};
    GrB_Index C[9] = {0, 1, 0, 1, 2, 2, 4, 3, 4};
    bool values[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    GrB_Matrix A = NULL;
    OK(GrB_Matrix_new(&A, GrB_BOOL, 5, 5));
    OK(GrB_Matrix_build_BOOL(A, R, C, values, 9, GrB_FIRST_BOOL));
    GxB_Matrix_fprint(A, "A", GxB_COMPLETE, stdout);
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    OK(LAGraph_MaximumMatching(&mateC, G, msg));
    printf("msg: %s\n", msg);

    LAGraph_Finalize(msg);
}

TEST_LIST =
    {
        {"Dummy", test_MCM}, // just one test in this example
        {NULL, NULL}};