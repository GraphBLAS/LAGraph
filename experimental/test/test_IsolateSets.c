#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include "LG_Xtest.h"

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G;
GrB_Matrix A;

#define dbg(x) GxB_print(x, 5)
#define err(x, info)                                    \
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) \
    {                                                   \
        char **err;                                     \
        GrB_error(err, x);                              \
        printf("\ninfo: %d error: %s\n", info, err);    \
    }
#define LEN 512
char filename[LEN + 1];
typedef struct
{
    const char *matrix_file; // Adjeancy matrix or graph
} matrix_info;
const matrix_info files[] = {
    {"empty.mtx"},
    {"comm0.mtx"},
    {"karate.mtx"},
    {""}};

void test_IsolateSets(void)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10_2
    LAGraph_Init(msg);
    printf("\n");
    // GrB_set(GrB_GLOBAL, 0, GxB_BURBLE);
    GrB_Info info;
    for (int k = 0;; k++)
    {
        if (strlen(files[k].matrix_file) == 0)
            break;
        snprintf(filename, LEN, LG_DATA_DIR "%s", files[k].matrix_file);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        OK(fclose(f));
        // GxB_print(A,5);
        GrB_Matrix MIset = NULL;


        double tsimple = LAGraph_WallClockTime();
        OK (LAGraph_IsolateSets(&MIset, A, 1231245, msg));
        tsimple = LAGraph_WallClockTime() - tsimple;
        printf(" time: %f\n", tsimple);

        if (MIset == NULL)
        {
            TEST_CHECK(true);
        }
        else
        {
            GrB_Index nrows;
            GrB_Matrix_nrows(&nrows, MIset);
            printf("%lu\n", nrows);
            GrB_Index nvals;
            GrB_Matrix_nvals(&nvals, MIset);
            printf("%lu\n", nvals);
            TEST_CHECK(nrows <= nvals);
        }
        // GxB_print(MIset,5);
        OK(LAGraph_Delete(&G, msg));
        GrB_Matrix_free(&A);
        GrB_Matrix_free(&MIset);
    }
    LAGraph_Finalize(msg);
#endif
};

TEST_LIST = {
    {"IsolateSets", test_IsolateSets},
    {NULL, NULL}};
