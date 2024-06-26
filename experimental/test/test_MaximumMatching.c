#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include "LG_internal.h"

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

#define LEN 512
char filename[LEN + 1];

void test_MCM()
{
    LAGraph_Init(msg);

    OK(LG_SET_BURBLE(1));

    GrB_Matrix A = NULL;
    snprintf(filename, LEN, LG_DATA_DIR "%s", "west0067.mtx");
    FILE *f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    GrB_Index nrows = 0, ncols = 0, nvals = 0;
    OK(GrB_Matrix_nrows(&nrows, A));
    OK(GrB_Matrix_ncols(&ncols, A));
    OK(GrB_Matrix_nvals(&nvals, A));

    // make A a bool matrix and iso-valued
    GrB_Index *I, *J, *X;
    double *dummy;
    OK(LAGraph_Malloc((void **)&I, nvals, sizeof(GrB_Index), msg));
    OK(LAGraph_Malloc((void **)&J, nvals, sizeof(GrB_Index), msg));
    OK(LAGraph_Malloc((void **)&dummy, nvals, sizeof(double), msg));

    GrB_Index Ibytes = 0, Jbytes = 0, Xbytes = 0;
    bool jumbled = 1;
    bool iso_value[nvals];
    for (uint64_t i = 0; i < nvals; i++)
        iso_value[i] = 1;
    OK(GrB_Matrix_extractTuples_FP64(I, J, dummy, &nvals, A));
    TEST_CHECK(I != NULL);
    OK(GrB_Matrix_new(&A, GrB_BOOL, nrows, ncols));
    OK(GrB_Matrix_build_BOOL(A, I, J, iso_value, nvals, GrB_FIRST_BOOL));

    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    TEST_CHECK(A == NULL); // A has been moved into G->A

    GrB_Vector mateC = NULL;
    OK(GrB_Vector_new(&mateC, GrB_UINT64, ncols));

    OK(LAGraph_MaximumMatching(&mateC, G, msg));
    printf("msg: %s\n", msg);

    GrB_Index nmatched = 0;

    GrB_Vector mateR = NULL;
    OK(GrB_Vector_new(&mateR, GrB_UINT64, nrows));

    // invert to check for dups
    OK(GxB_Vector_unpack_CSC(mateC, (GrB_Index **)&J, (void **)&X, &Jbytes, &Xbytes, NULL, &nmatched, &jumbled, NULL));
    OK(GrB_Vector_build_UINT64(mateR, X, J, nmatched, GrB_FIRST_UINT64));
    GrB_Index nmateR = 0;
    OK(GrB_Vector_nvals(&nmateR, mateR));
    // if nvals of mateC and mateR don't match, then there's at least one row that is used in at least one matching
    TEST_CHECK(nmatched == nmateR);
    OK(LAGraph_Free((void **)&mateR, msg));

    // pack matched values in a matrix
    GrB_Matrix M = NULL;
    A = G->A;
    bool val[nmatched];
    for (uint64_t i = 0; i < nmatched; i++)
        val[i] = 1;
    OK(GrB_Matrix_new(&M, GrB_BOOL, nrows, ncols));
    OK(GrB_Matrix_build_BOOL(M, X, J, val, nmatched, NULL));
    // mask with matrix A to check if all edges are present in A
    OK(GrB_Matrix_assign(M, M, NULL, A, GrB_ALL, nrows, GrB_ALL, ncols, GrB_DESC_S));
    GrB_Index nvalsM = 0;
    OK(GrB_Matrix_nvals(&nvalsM, M));
    // if values have been eliminated then edges do not exist in A
    TEST_CHECK(nvalsM == nmatched);

    // sprank must be equal to nvals of mateC (nmatched)

    LAGraph_Finalize(msg);
}

TEST_LIST =
    {
        {"MaximumMatching", test_MCM}, // just one test in this example
        {NULL, NULL}};