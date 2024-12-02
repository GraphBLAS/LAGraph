#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include "LG_Xtest"
#include <stdio.h>


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
GrB_Index S, T;
#define LEN 512
#define NTESTS 1
char filename[LEN + 1];

const char* files[NTESTS] = {}; // add matrix files

//maybe add arrays for T, S, and F??

void test_MaxFlow(void) {
  LAGraph_Init(msg);
  OK(LG_SET_BURBLE(1));
  
  for(uint8_t jit = 0; jit < 2; jit++){
    uint8_t JIT_flag = jit * 4;
    OK(GxB_Gloabl_Option_set(GxB_JIT_C_CONTROL, JIT_flag));
    for(uint8_t test = 0; test < NTESTS; test++){
      GrB_Matrix A;
      GrB_Index S, T;
      snprintf(filename, LEN, LG_DATA_DIR "%s", filenames[test]);
      FILE* f = fopen(filename, "r");
      TEST_CHECK(f != NULL);
      OK(LAGraph_MMRead(&A, f, msg));
      OK(fclose(f));
      GrB_Index nrows = 0, ncols = 0, nvals = 0;
      OK(GrB_Matrix_nrows(&nrows, A));
      OK(GrB_Matrix_ncols(&ncols, A));
      OK(GrB_Matrix_nvals(&nvlas, A));

      GrB_Index *I, *J;
      float * vals;
      OK(LAGraph_Malloc((void**)&I, nvals, sizeof(GrB_Index), msg));
      OK(LAGraph_Malloc((void**)&J, nvals, sizeof(GrB_Index), msg));
      OK(LAGraph_Malloc((void**)&vals, nvals, sizeof(GrB_FP32), msg));

      OK(GrB_Matrix_extractTuples(I, J, dummy, &nvals, A));
      TEST_CHECK(I != NULL);
      OK(GrB_Matrix_new(&A, GrB_FP32, nrows, ncols));
      OK(GrB_Matrix_build(A, I, J, vals, nvals, GxB_FIRST_FP32));

      OK(LAGraph_free((void**)&I, msg));
      OK(LAGraph_free((void**)&J, msg));
      OK(LAGraph_free((void**)&vals, msg));
    }
  }
}
