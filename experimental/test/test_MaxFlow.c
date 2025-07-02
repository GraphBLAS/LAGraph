//------------------------------------------------------------------------------
// experimental/test/test_MaxFlow: tests for LAGr_MaxFlow
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Darin Peries and Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <stdio.h>
#include <LG_Xtest.h>
#include "LG_internal.h"


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
#define LEN 512
#define NTESTS 7
char filename[LEN + 1];

typedef struct{
  char* filename;
  GrB_Index S;
  GrB_Index T;
  double F;
}test_info;

test_info tests[] = {
  {"wiki.mtx", 0, 5, 4},
  {"matrix_random_flow.mtx", 0,9, 22},
  {"rand.mtx", 0, 19, 37}, 
  {"mcl.mtx", 0, 9, 0}, 
  {"cycle_flow.mtx", 0, 89, 1},
  {"random_weighted_general2.mtx", 0, 299, 11098623877},
  {"random_weighted_general1.mtx", 0, 499, 6264009335}
};

//399 11098623877 alt sink and src for test 6

void test_MaxFlow(void) {
#if LG_SUITESPARSE_GRAPHBLAS_V10
  LAGraph_Init(msg);
//OK(LG_SET_BURBLE(1));
  OK(LG_SET_BURBLE(0));
  OK(GxB_Global_Option_set(GxB_JIT_C_CONTROL, 4));
  for(uint8_t test = 0; test < NTESTS; test++){
    GrB_Matrix A=NULL;
    TEST_CASE(tests[test].filename);
    snprintf(filename, LEN, LG_DATA_DIR "%s", tests[test].filename);
    FILE* f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    OK(LAGraph_Cached_AT(G, msg));
    OK(LAGraph_Cached_EMin(G, msg));

    //begin test
    double flow = 0;
    OK(LAGr_MaxFlow(&flow, NULL, G, tests[test].S, tests[test].T, msg));
    printf("%s\n", msg);
    TEST_CHECK(flow == tests[test].F);
    printf("flow is: %lf\n", flow);

    //free work
    OK(LAGraph_Delete(&G, msg));
  }
  LAGraph_Finalize(msg);
#endif
}

void test_MaxFlowMtx(void) {
  LAGraph_Init(msg);
#if LG_SUITESPARSE_GRAPHBLAS_V10
//OK(LG_SET_BURBLE(1));
  OK(LG_SET_BURBLE(0));
  OK(GxB_Global_Option_set(GxB_JIT_C_CONTROL, 4));
  for(uint8_t test = 0; test < NTESTS; test++){
    GrB_Matrix A=NULL;
    TEST_CASE(tests[test].filename);
    snprintf(filename, LEN, LG_DATA_DIR "%s", tests[test].filename);
    FILE* f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));

    //create flow mtx
    GrB_Matrix flow_mtx=NULL;
    GrB_Index n;
    OK(GrB_Matrix_nrows(&n, A));
    OK(GrB_Matrix_new(&flow_mtx, GrB_FP64, n, n));
    
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    OK(LAGraph_Cached_AT(G, msg));
    OK(LAGraph_Cached_EMin(G, msg));

    //begin test
    double flow = 0;
    OK(LAGr_MaxFlow(&flow, &flow_mtx, G, tests[test].S, tests[test].T, msg));
    int status = LG_check_flow(&flow_mtx, msg);
    printf("%d", status);
    printf("%s\n", msg);
    TEST_CHECK(flow == tests[test].F);
    printf("flow is: %lf\n", flow);

    //free work
    GrB_free(&flow_mtx);
    OK(LAGraph_Delete(&G, msg));
  }
#endif
  LAGraph_Finalize(msg);
}


TEST_LIST = {{"MaxFlow", test_MaxFlow}, {"MaxFlowMtx", test_MaxFlowMtx}, {NULL, NULL}};
