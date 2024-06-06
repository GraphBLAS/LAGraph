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

    GrB_Matrix A = NULL;
    OK(GrB_Matrix_new(&A, GrB_BOOL, 5, 5));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    OK(LAGraph_MaximumMatching(&mateC, G, msg));
    printf("msg: %s\n", msg);

    LAGraph_Finalize(msg);
}

TEST_LIST =
    {
        {"Dummy", test_MCM}, // just one test in this example
        {NULL, NULL}};