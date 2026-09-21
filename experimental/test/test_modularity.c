
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
    {"comm0.mtx", "comm0_S.mtx", 1, -0.17347},
    {"comm0.mtx", "comm0_Sa.mtx", 1, 0.35714},
    // {"com-Amazon.mtx", "comm0_Sa.mtx",1, -1},
    {"", "", -1, -1}};

typedef struct
{
    const char *matrix_file;
    const char *cluster_file;
    const double gamma;
} large_matrix_info;

const large_matrix_info large_files[] = {
    {"jagmesh7.mtx", "jagmesh7_cluster.mtx", 1},
    {"bcsstk13.mtx", "bcsstk13_cluster.mtx", 1},
    {"west0067.mtx", "west0067_cluster.mtx", 1},
    {"", "", -1}};

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
        OK(LAGr_AdjModularity(&Q, gamma, A, S, msg));
        Q = floor(100000*Q)/100000;
        TEST_CHECK(Q == files[k].mod);
        printf("Q:%.15g\n", Q);
        OK(GrB_free(&A));
        OK(GrB_free(&S));
    }

    for (int k = 0;; k++)
    {
        if (strlen(large_files[k].matrix_file) == 0)
            break;

        snprintf(filename, LEN, LG_DATA_DIR "%s", large_files[k].matrix_file);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        OK(fclose(f));

        GrB_Vector c = NULL;
        snprintf(filename2, LEN, LG_DATA_DIR "%s", large_files[k].cluster_file);
        FILE *t = fopen(filename2, "r");
        TEST_CHECK(t != NULL);
        OK(LAGraph_MMRead((GrB_Matrix *)&c, t, msg));
        OK(fclose(t));

        GrB_Index n = 0, csize = 0;
        OK(GrB_Matrix_nrows(&n, A));
        OK(GrB_Vector_size(&csize, c));
        TEST_CHECK(n == csize);

        GrB_Index cnvals = 0;
        OK(GrB_Vector_nvals(&cnvals, c));

        GrB_Index *I = NULL;
        int64_t *X = NULL;
        int64_t *labels = NULL;
        int64_t *unique_labels = NULL;
        GrB_Index *comm_ids = NULL;

        OK(LAGraph_Malloc((void **)&I, cnvals, sizeof(GrB_Index), msg));
        OK(LAGraph_Malloc((void **)&X, cnvals, sizeof(int64_t), msg));
        OK(LAGraph_Calloc((void **)&labels, n, sizeof(int64_t), msg));
        OK(LAGraph_Malloc((void **)&unique_labels, n, sizeof(int64_t), msg));
        OK(LAGraph_Malloc((void **)&comm_ids, n, sizeof(GrB_Index), msg));

        OK(GrB_Vector_extractTuples_INT64(I, X, &cnvals, c));
        for (GrB_Index p = 0; p < cnvals; p++)
        {
            labels[I[p]] = X[p];
        }

        GrB_Index ncomms = 0;
        for (GrB_Index i = 0; i < n; i++)
        {
            int64_t label = labels[i];
            GrB_Index cidx = 0;
            bool found = false;
            for (GrB_Index j = 0; j < ncomms; j++)
            {
                if (unique_labels[j] == label)
                {
                    cidx = j;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                unique_labels[ncomms] = label;
                cidx = ncomms;
                ncomms++;
            }
            comm_ids[i] = cidx;
        }

        OK(GrB_Matrix_new(&S, GrB_BOOL, n, n));
        for (GrB_Index i = 0; i < n; i++)
        {
            OK(GrB_Matrix_setElement_BOOL(S, true, i, comm_ids[i]));
        }

        double Q;
        OK(LAGr_AdjModularity(&Q, large_files[k].gamma, A, S, msg));
        TEST_CHECK(isfinite(Q));

        OK(LAGraph_Free((void **)&I, msg));
        OK(LAGraph_Free((void **)&X, msg));
        OK(LAGraph_Free((void **)&labels, msg));
        OK(LAGraph_Free((void **)&unique_labels, msg));
        OK(LAGraph_Free((void **)&comm_ids, msg));
        OK(GrB_free(&c));
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
