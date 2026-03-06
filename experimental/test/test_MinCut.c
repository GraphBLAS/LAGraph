#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
#define LEN 512
#define NTESTS 7
char filename[LEN + 1];


typedef struct {
  char* filename;
  GrB_Index s;
  GrB_Index t;
  LAGraph_Kind kind;
} test_case;

test_case tests[] = {
  {"wiki.mtx", 0, 5, LAGraph_ADJACENCY_DIRECTED},
  {"matrix_random_flow.mtx", 0,9, LAGraph_ADJACENCY_DIRECTED},
  {"rand.mtx", 0, 19, LAGraph_ADJACENCY_DIRECTED},
  {"mcl.mtx", 0, 9, LAGraph_ADJACENCY_DIRECTED},
  {"cycle_flow.mtx", 0, 89,  LAGraph_ADJACENCY_DIRECTED},
  {"random_weighted_general2.mtx", 0, 299, LAGraph_ADJACENCY_UNDIRECTED},
  {"random_weighted_general1.mtx", 0, 499, LAGraph_ADJACENCY_UNDIRECTED}
};


void test_MinCut() {

  LAGraph_Init(msg);
#if LG_SUITESPARSE_GRAPHBLAS_V10
  OK(LG_SET_BURBLE(0));
  for(uint8_t test = 0; test < NTESTS; test++){
    GrB_Matrix A=NULL, R=NULL, cut_set=NULL;
    GrB_Vector S=NULL, S_bar=NULL;
    double min_cut = 0;
    GrB_Index n = 0;
    printf ("\nMatrix: %s\n", tests[test].filename);
    TEST_CASE(tests[test].filename);
    
    snprintf(filename, LEN, LG_DATA_DIR "%s", tests[test].filename);
    FILE* f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    
    OK(LAGraph_MMRead(&A, f, msg));
    
    OK(fclose(f));
    LAGraph_Kind kind = tests [test].kind ;
    OK(LAGraph_New(&G, &A, kind, msg));
    if (kind == LAGraph_ADJACENCY_DIRECTED)
    {
        OK(LAGraph_Cached_AT(G, msg));
    }

    OK(LAGraph_Cached_EMin(G, msg));

    // test with JIT
    OK(GxB_Global_Option_set(GxB_JIT_C_CONTROL, GxB_JIT_ON));
    double flow = 0;
    OK(LAGr_MaxFlow(&flow, NULL, &R, G, tests[test].s, tests[test].t, msg));
    printf("%s\n", msg);
    printf("flow is: %lf\n", flow);

    OK(LAGraph_MinCut(&S, &S_bar, &cut_set, G, R, tests[test].s, tests[test].t, msg));
    printf("%s\n", msg);

    OK(GrB_reduce(&min_cut, NULL, GrB_PLUS_MONOID_FP64, cut_set, NULL));

    TEST_CHECK(flow == min_cut);
    printf("The min cut: %lf\n", min_cut);

    GrB_free(&A);
    GrB_free(&R);
    GrB_free(&S);
    GrB_free(&S_bar);
    GrB_free(&cut_set);
    LAGraph_Delete(&G, msg);
  }

#endif
  

  LAGraph_Finalize(msg);
  
}

TEST_LIST = {{"MinCut", test_MinCut}, {NULL, NULL}};
