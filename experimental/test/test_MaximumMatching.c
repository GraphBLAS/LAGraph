#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include "LG_internal.h"


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

void test_MCM()
{
    LAGraph_Init(msg);

    // OK(LG_SET_BURBLE(1));

    GrB_Index R[9] = {0, 0, 1, 2, 2, 3, 3, 4, 4};
    GrB_Index C[9] = {0, 1, 0, 1, 2, 2, 4, 3, 4};
    bool values[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    GrB_Matrix A = NULL;
    GrB_Index nrows = 5, ncols = 5;
    OK(GrB_Matrix_new(&A, GrB_BOOL, nrows, ncols));
    OK(GrB_Matrix_build_BOOL(A, R, C, values, 9, GrB_FIRST_BOOL)); // change to pack
    // GxB_Matrix_fprint(A, "A", GxB_COMPLETE, stdout);
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)); // A is set to NULL afterwards

    GrB_Vector mateC = NULL;
    OK(GrB_Vector_new(&mateC, GrB_UINT64, ncols));
    GrB_Index Ilist[2] = {2, 3};
    uint64_t Vlist[2] = {2, 4};
    OK(GrB_Vector_build_UINT64(mateC, Ilist, Vlist, 2, NULL));

    OK(LAGraph_MaximumMatching(&mateC, G, msg));
    printf("msg: %s\n", msg);

    GrB_Index *J = NULL, *X = NULL;
    GrB_Index Jbytes = 0, Xbytes = 0, nmatched = 0;
    bool jumbled = 1;

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
    bool x[1] = {1};
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
        {"Dummy", test_MCM}, // just one test in this example
        {NULL, NULL}};