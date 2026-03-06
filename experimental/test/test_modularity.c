
#include <acutest.h>
#include <stdio.h>

#include "LG_Xtest.h"
#include <LAGraphX.h>
#include <LAGraph_test.h>
char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
GrB_Matrix S = NULL;
#define LEN 512
char filename[LEN + 1];
char filename2[LEN + 1];
typedef struct
{
    const char *matrix_file; // Adjeancy matrix or graph
    const char *comm_matrix; // comm matrix define as pattern as general
    const double gamma;      // resolution of communities 1 for now
    const double mod;        // expected modularity from Q = 1/2m (S^TBS)
} matrix_info;

const matrix_info files[] = {   
    {"empty.mtx","empty.mtx",1,0},
    {"comm0.mtx", "comm0_S.mtx", 1, -0.17347},
    {"comm0.mtx", "comm0_Sa.mtx", 1, 0.35714},
    // {"com-Amazon.mtx", "comm0_Sa.mtx",1, -1},
    {"", "", -1, -1}};
void test_modularity(void)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10_2
    LAGraph_Init(msg);
    for (int k = 0;; k++)
    {
        if (strlen(files[k].matrix_file) == 0)
            break;
        snprintf(filename, LEN, LG_DATA_DIR "%s", files[k].matrix_file);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        OK(fclose(f));
        printf("\nInput of Matrix:\n");
        GxB_print(A, 3);

        snprintf(filename2, LEN, LG_DATA_DIR "%s", files[k].comm_matrix);
        FILE *t = fopen(filename2, "r");
        TEST_CHECK(t != NULL);
        OK(LAGraph_MMRead(&S, t, msg));
        OK(fclose(t));
        printf("\nInput of Matrix S:\n");
        GxB_print(S, 3);

        double gamma = files[k].gamma;
        double Q;
        OK(LAGr_ModularityMatrix(&Q, gamma, A, S, msg));
        Q = floor(100000*Q)/100000;
        TEST_CHECK(Q == files[k].mod);
        printf("Q:%.15g\n", Q);
        OK(GrB_free(&A));
        OK(GrB_free(&S));
    }

    OK(GrB_free(&A));
    OK(GrB_free(&S));

    OK(LAGraph_Delete(&G, msg));

    LAGraph_Finalize(msg);
    #endif
}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
    {
        {"Modularity", test_modularity}, // just one test in this example
        {NULL, NULL}};
