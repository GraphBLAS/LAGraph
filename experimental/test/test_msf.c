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
    double ans_w;
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
    { 1, "A.mtx", 6, A_mtx_i, A_mtx_j, NAN},
    { 1, "jagmesh7.mtx", 1137, NULL, NULL, NAN},
    { 0, "west0067.mtx", 66, NULL, NULL, -63.9103636}, // unsymmetric
    { 1, "bcsstk13.mtx", 2002, NULL, NULL, -27812381075940.4},
    { 0, "matrix_int8.mtx", 6, mtx8_i, mtx8_j, -120.0},
    { 0, "matrix_uint8.mtx", 6, mtx8u_i, mtx8u_j, 8.0},
    { 1, "karate.mtx", 33, NULL, NULL, NAN},
    { 1, "ldbc-cdlp-undirected-example.mtx", 7, NULL, NULL, NAN},
    { 1, "ldbc-undirected-example-bool.mtx", 8, NULL, NULL, NAN},
    { 1, "ldbc-undirected-example-unweighted.mtx", 8, NULL, NULL, NAN},
    { 1, "ldbc-undirected-example.mtx", 8, NULL, NULL, NAN},
    { 1, "ldbc-wcc-example.mtx", 9, NULL, NULL, NAN},
    { 0, "" },
} ;

//****************************************************************************
void test_msf (void)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LAGraph_Init (msg) ;
    bool burble = false ;
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

        OK (LAGraph_Matrix_Print (A, LAGraph_SHORT, stdout, msg)) ;
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
            if (jit) printf ("\nJIT is enabled\n") ; else printf ("\nJIT is disabled\n") ;

            // connected components
            GrB_Vector cc0 = NULL, cc1 = NULL, cc2 = NULL;

            OK (GxB_Global_Option_set (GxB_JIT_C_CONTROL,
                jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;
            // compute the min spanning forest
            C = NULL ;
            OK (LG_SET_BURBLE (burble)) ;
            int result = LAGraph_msf (&C, &cc0, A, sanitize, msg) ;
            OK (LG_SET_BURBLE (false)) ;

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
            OK(LAGraph_New(&G, &S, LAGraph_ADJACENCY_UNDIRECTED, msg));
 
            //Check that the graph has all the same ccs.
            OK (LAGr_ConnectedComponents(&cc1, G, msg)) ;
            bool ok = false ;
            OK (GrB_Vector_new(&cc2, GrB_UINT64, n)) ;
            // cc1 and cc0 should have the same structure as cc2. 
            // make their values equal and then compare them.
            // msf does not guarentee that the lower node is used as componentId
            OK (GxB_Vector_extract_Vector(cc2, NULL, NULL, cc0, cc1, NULL)) ;
            OK (LAGraph_Vector_IsEqual(&ok, cc2, cc0, msg)) ;

//          if(!ok)
//          {
//              GxB_print(cc2, GxB_SHORT);
//              GxB_print(cc0, GxB_SHORT);
//          }
            TEST_ASSERT(ok) ;
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
            double tot_weight, ans_w = files[k].ans_w;
            OK (GrB_Matrix_reduce_FP64(&tot_weight, NULL, GrB_PLUS_MONOID_FP64, C, NULL)) ;
            TEST_CHECK (isnan(files[k].ans_w) || 
                    fabs(tot_weight - ans_w) <= 1E-10 * fabs(ans_w)) ;           
            OK (LAGraph_Matrix_Print (C, pr, stdout, msg)) ;
            OK (LAGraph_Delete (&G, msg)) ;
            OK (GrB_free (&cc0)) ;
            OK (GrB_free (&cc1)) ;
            OK (GrB_free (&cc2)) ;
            OK (GrB_free (&C)) ;

            printf ("JIT test is done\n") ;
        }
        OK (GrB_free(&Ans)) ;
        OK (GrB_free (&A)) ;
    }
    GrB_free(&zeroB);
    LAGraph_Finalize (msg) ;
    #endif
}

//------------------------------------------------------------------------------
// infinity test
//------------------------------------------------------------------------------

void test_inf_msf (void)
{
    LAGraph_Init (msg) ;
    bool burble = false ;
    GrB_Scalar zeroB = NULL;
    GrB_Scalar_new(&zeroB, GrB_BOOL);
    GrB_Scalar_setElement_BOOL(zeroB, false);

    // load the matrix as A
    const char *aname = "bcsstk13.mtx" ;
    bool symmetric = 1 ;
    printf ("\n================================== %s:\n", aname) ;
    TEST_CASE (aname) ;
    snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    fclose (f) ;

    GrB_Index n = 0;
    OK (GrB_Matrix_nrows (&n, A)) ;
    OK (GrB_Matrix_new(&S, GrB_BOOL, n, n)) ;

    // Make A be iso infinity and S iso true.
    OK (GrB_Matrix_assign_FP64(
        A, A, NULL, INFINITY, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
    OK (GrB_Matrix_assign_BOOL(
        S, A, NULL, (bool) true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;

    // compute the min spanning forest
    S_C = C = NULL ;
    OK (LG_SET_BURBLE (burble)) ;
    int result = LAGraph_msf (&C, NULL, A, false, msg) ;
    OK (LG_SET_BURBLE (false)) ;
    printf ("result: %d\n", result) ;
    OK(result);

    OK (LG_SET_BURBLE (burble)) ;
    result = LAGraph_msf (&S_C, NULL, S, false, msg) ;
    OK (LG_SET_BURBLE (false)) ;
    printf ("result: %d\n", result) ;
    OK(result);

    bool ok = false ;
    // check structure is equal.
    OK (LAGraph_Matrix_IsEqualOp(&ok, C, S_C, GrB_ONEB_BOOL, msg));
    TEST_CHECK(ok);
    OK (GrB_free (&S)) ;
    OK (GrB_free (&S_C)) ;
    OK (GrB_free (&C)) ;
    OK (GrB_free (&A)) ;
    GrB_free(&zeroB);
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors (void)
{
    LAGraph_Init (msg) ;

    #if LG_SUITESPARSE_GRAPHBLAS_V10
    // C and A are NULL
    int result = LAGraph_msf (NULL, NULL, NULL, true, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // A must be square
    OK (GrB_Matrix_new (&A, GrB_UINT64, 3, 4)) ;
    result = LAGraph_msf (&C, NULL, A, true, msg) ;
    TEST_CHECK (result == GrB_DIMENSION_MISMATCH) ;
    OK (GrB_free (&A)) ;

    // A must real
    OK (GrB_Matrix_new (&A, GxB_FC32, 4, 4)) ;
    result = LAGraph_msf (&C, NULL, A, true, msg) ;
    TEST_CHECK (result == GrB_DOMAIN_MISMATCH) ;
    
    #else 
    // Not implemented
    OK (GrB_Matrix_new (&A, GrB_BOOL, 4, 4)) ;
    int result = LAGraph_msf (&C, NULL, A, true, msg) ;
    TEST_CHECK (result == GrB_NOT_IMPLEMENTED) ;
    #endif

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

//****************************************************************************

TEST_LIST = {
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    {"msf", test_msf},
    {"inf_msf", test_inf_msf},
    #endif
    {"msf_errors", test_errors},
    {NULL, NULL}
};
