//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_msf.c: test cases for Min Spanning Forest
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// todo: write a simple msf method, as LG_check_msf, and compare its results
// with LAGraph_msf

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL, G_C = NULL;
GrB_Matrix A = NULL ;
GrB_Matrix S = NULL ;
GrB_Matrix S_C = NULL ;
GrB_Matrix C = NULL ;
GrB_Matrix Ans = NULL ;
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    bool symmetric ;
    const char *name ;
    uint64_t ans_n;
    const uint64_t *ans_i;
    const uint64_t *ans_j;
}
matrix_info ;
const uint64_t A_mtx_i [] = {1, 2, 3, 4, 5, 6};
const uint64_t A_mtx_j [] = {0, 0, 1, 1, 1, 0};
const uint64_t mtx8u_i [] = {1, 2, 3, 4, 5, 6};
const uint64_t mtx8u_j [] = {4, 3, 0, 6, 2, 3};
const uint64_t mtx8_i [] = {1, 2, 3, 4, 5, 6};
const uint64_t mtx8_j [] = {4, 3, 0, 6, 2, 3};
const matrix_info files [ ] =
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10 
    { 1, "A.mtx", 6, A_mtx_i, A_mtx_j},
    { 1, "jagmesh7.mtx", 1137, NULL, NULL},
    { 0, "west0067.mtx", 66, NULL, NULL}, // unsymmetric
    { 1, "bcsstk13.mtx", 2002, NULL, NULL}, // overflows an INT32
    { 0, "matrix_int8.mtx", 6, mtx8_i, mtx8_j},
    { 0, "matrix_uint8.mtx", 6, mtx8u_i, mtx8u_j},
    { 1, "karate.mtx", 33, NULL, NULL},
    { 1, "ldbc-cdlp-undirected-example.mtx", 7, NULL, NULL},
    { 1, "ldbc-undirected-example-bool.mtx", 8, NULL, NULL},
    { 1, "ldbc-undirected-example-unweighted.mtx", 8, NULL, NULL},
    { 1, "ldbc-undirected-example.mtx", 8, NULL, NULL},
    { 1, "ldbc-wcc-example.mtx", 9, NULL, NULL},
    #endif 
    { 0, "" },
} ;

//****************************************************************************
void test_msf (void)
{
    #if LAGRAPH_SUITESPARSE
    LAGraph_Init (msg) ;
    GrB_Scalar zeroB = NULL;
    GrB_Scalar_new(&zeroB, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(zeroB, false);
    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        bool symmetric = files [k].symmetric ;
        uint64_t branches = 0;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;

        // ensure A is uint64
        GrB_Index n = 0;
        OK (GrB_Matrix_nrows (&n, A)) ;

        // construct a directed graph G with adjacency matrix S
        TEST_CHECK (S == NULL) ;

        OK (LAGraph_Matrix_Print (A, GxB_SHORT, stdout, msg)) ;
        bool sanitize = (!symmetric) ;

        if (files[k].ans_i && files[k].ans_j)
        {
            OK (GrB_Matrix_new(&Ans, GrB_BOOL, n, n)) ;
            OK (GxB_Matrix_build_Scalar(
                Ans, files[k].ans_i, files[k].ans_j, zeroB, files[k].ans_n
            )) ;
        }
        
        for (int jit = 0 ; jit <= 1 ; jit++)
        {
            OK (GxB_Global_Option_set (GxB_JIT_C_CONTROL,
                jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;
            // compute the min spanning forest
            C = NULL ;
            // GxB_Global_Option_set(GxB_BURBLE, true);
            int result = LAGraph_msf (&C, A, sanitize, msg) ;
            // GxB_Global_Option_set(GxB_BURBLE, false);
            printf ("result: %d\n", result) ;
            OK(result);
            GrB_Matrix_nvals(&branches, C);
            TEST_CHECK(branches == files[k].ans_n);
            LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT ;
            
            OK (GrB_Matrix_new(&S, GrB_BOOL, n, n)) ;
            OK (GrB_Matrix_assign_BOOL(
                S, A, NULL, (bool) true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
            if(!symmetric)
            {
                OK (GrB_Matrix_eWiseAdd_BinaryOp(
                    S, NULL, NULL, GxB_ANY_BOOL, S, S, GrB_DESC_T1)) ;
            }
            OK (GrB_Matrix_new(&S_C, GrB_BOOL, n, n)) ;
            OK (GrB_Matrix_assign_BOOL(
                S_C, C, NULL, (bool) true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
            OK (GrB_Matrix_eWiseAdd_BinaryOp(
                S_C, NULL, NULL, GxB_ANY_BOOL, S_C, S_C, GrB_DESC_T1)) ;
            OK(LAGraph_New(&G, &S, LAGraph_ADJACENCY_UNDIRECTED, msg));
            OK(LAGraph_New(&G_C, &S_C, LAGraph_ADJACENCY_UNDIRECTED, msg));

            
            //Check that the graph has all the same ccs.
            GrB_Vector cc0 = NULL, cc1 = NULL;
            OK (LAGr_ConnectedComponents(&cc0, G, msg));
            OK (LAGr_ConnectedComponents(&cc1, G_C, msg));
            bool ok = false ;
            OK (LAGraph_Vector_IsEqual(&ok, cc0, cc1, msg));
            TEST_CHECK(ok);
            // check result C for A.mtx
            if (files[k].ans_i && files[k].ans_j)
            {
                OK (GrB_Matrix_eWiseMult_BinaryOp(
                    Ans, NULL, GxB_LOR_BOOL, GrB_ONEB_BOOL, Ans, C, NULL)) ;
                OK (GrB_Matrix_eWiseMult_BinaryOp(
                    Ans, NULL, GxB_LOR_BOOL, GrB_ONEB_BOOL, Ans, C, GrB_DESC_T1
                )) ;
                OK (GrB_Matrix_reduce_BOOL(
                    &ok, NULL, GrB_LAND_MONOID_BOOL, Ans, NULL));
                TEST_CHECK (ok) ;
            }

            printf ("\nmsf:\n") ;
            OK (LAGraph_Matrix_Print (C, pr, stdout, msg)) ;
            OK (LAGraph_Delete (&G, msg)) ;
            OK (LAGraph_Delete (&G_C, msg)) ;
            OK (GrB_free (&cc0)) ;
            OK (GrB_free (&cc1)) ;
            OK (GrB_free (&C)) ;
        }
        OK (GrB_free(&Ans)) ;
        OK (GrB_free (&A)) ;
    }
    GrB_free(&zeroB);
    LAGraph_Finalize (msg) ;
    #endif
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors (void)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LAGraph_Init (msg) ;

    // C and A are NULL
    int result = LAGraph_msf (NULL, NULL, true, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // A must be square
    OK (GrB_Matrix_new (&A, GrB_UINT64, 3, 4)) ;
    result = LAGraph_msf (&C, A, true, msg) ;
    TEST_CHECK (result == GrB_DIMENSION_MISMATCH) ;
    OK (GrB_free (&A)) ;

    // A must be square
    OK (GrB_Matrix_new (&A, GxB_FC32, 4, 4)) ;
    result = LAGraph_msf (&C, A, true, msg) ;
    TEST_CHECK (result == GrB_DOMAIN_MISMATCH) ;

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
    #endif
}

//****************************************************************************

TEST_LIST = {
    {"msf", test_msf},
    {"msf_errors", test_errors},
    {NULL, NULL}
};
