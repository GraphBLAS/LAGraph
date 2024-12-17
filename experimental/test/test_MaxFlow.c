#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <stdio.h>

#include "LG_Xtest.h"
#include "LG_internal.h"


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
#define LEN 512
#define NTESTS 2
char filename[LEN + 1];

typedef struct{
  char* filename;
  GrB_Index S;
  GrB_Index T;
  float F;
}test_info;

test_info tests[] = {
  {"wiki.mtx", 0, 5, 4},
  {"matrix_random_flow.mtx", 0,9, 22}
};


void test_MaxFlow(void) {
  LAGraph_Init(msg);
  OK(LG_SET_BURBLE(1));
  OK(GxB_Global_Option_set(GxB_JIT_C_CONTROL, 4));
  for(uint8_t test = 0; test < NTESTS; test++){
    GrB_Matrix A;
    TEST_CASE(tests[test].filename);
    snprintf(filename, LEN, LG_DATA_DIR "%s", tests[test].filename);
    FILE* f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    //OK(GrB_Matrix_new(&A, GrB_FP32, nrows, ncols));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));

    //begin test
    float flow = 0;
    OK(LAGraph_MaxFlow(G, tests[test].S, tests[test].T, &flow, msg));
    printf("%s\n", msg);
    TEST_CHECK(flow == tests[test].F);
    printf("flow is: %f\n", flow);

    //free work
    GrB_free(&A);
    OK(LAGraph_Delete(&G, msg));
  }
  LAGraph_Finalize(msg);
}

TEST_LIST = {{"MaxFlow", test_MaxFlow}, {NULL, NULL}};
