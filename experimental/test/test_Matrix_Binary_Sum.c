//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_Matrix_Binary_Sum.c:  test LAGraph_Matrix_Binary_Sum
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Michel Pelletier.

//------------------------------------------------------------------------------

#include "LAGraph_test.h"
#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL, B = NULL, C = NULL, Expected = NULL ;

//------------------------------------------------------------------------------
// setup and teardown
//------------------------------------------------------------------------------

void setup (void)
{
    OK (LAGraph_Init (msg)) ;
}

void teardown (void)
{
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// make_AB: construct two 4x4 FP64 matrices with overlapping and distinct
// (i,j) coordinates
//------------------------------------------------------------------------------

static void make_AB (void)
{
    // A has entries at (0,0), (1,1), (2,2), (0,3)
    GrB_Index Ai [ ] = { 0, 1, 2, 0 } ;
    GrB_Index Aj [ ] = { 0, 1, 2, 3 } ;
    double    Ax [ ] = { 1, 2, 3, 4 } ;
    OK (GrB_Matrix_new (&A, GrB_FP64, 4, 4)) ;
    OK (GrB_Matrix_build_FP64 (A, Ai, Aj, Ax, 4, NULL)) ;

    // B has entries at (0,0), (1,1), (3,3), (2,0)
    //   shares (0,0) and (1,1) with A, distinct elsewhere
    GrB_Index Bi [ ] = { 0, 1, 3, 2 } ;
    GrB_Index Bj [ ] = { 0, 1, 3, 0 } ;
    double    Bx [ ] = { 10, 20, 30, 40 } ;
    OK (GrB_Matrix_new (&B, GrB_FP64, 4, 4)) ;
    OK (GrB_Matrix_build_FP64 (B, Bi, Bj, Bx, 4, NULL)) ;
}

//------------------------------------------------------------------------------
// test_Matrix_Binary_Sum: basic correctness
//------------------------------------------------------------------------------

void test_Matrix_Binary_Sum (void)
{
    setup ( ) ;

    make_AB ( ) ;

    // Expected = A + B (element-wise add, set union)
    OK (GrB_Matrix_new (&Expected, GrB_FP64, 4, 4)) ;
    OK (GrB_eWiseAdd (Expected, NULL, NULL, GrB_PLUS_FP64, A, B, NULL)) ;

    // C = binary sum of {A, B} using GrB_PLUS_FP64 to combine entries
    GrB_Matrix Mats [2] = { A, B } ;
    OK (LAGraph_Matrix_Binary_Sum (&C, Mats, 2, GrB_PLUS_FP64, msg)) ;

    bool ok ;
    OK (LAGraph_Matrix_IsEqual (&ok, C, Expected, msg)) ;
    TEST_CHECK (ok) ;
    TEST_MSG ("binary sum of {A,B} did not equal A+B") ;
    OK (GrB_free (&C)) ;

    // single-matrix array: result must be a copy of A
    GrB_Matrix One [1] = { A } ;
    OK (LAGraph_Matrix_Binary_Sum (&C, One, 1, GrB_PLUS_FP64, msg)) ;
    OK (LAGraph_Matrix_IsEqual (&ok, C, A, msg)) ;
    TEST_CHECK (ok) ;
    TEST_MSG ("binary sum of {A} did not equal A") ;
    OK (GrB_free (&C)) ;

    // an empty matrix in the array contributes nothing
    GrB_Matrix Empty = NULL ;
    OK (GrB_Matrix_new (&Empty, GrB_FP64, 4, 4)) ;
    GrB_Matrix WithEmpty [3] = { A, Empty, B } ;
    OK (LAGraph_Matrix_Binary_Sum (&C, WithEmpty, 3, GrB_PLUS_FP64, msg)) ;
    OK (LAGraph_Matrix_IsEqual (&ok, C, Expected, msg)) ;
    TEST_CHECK (ok) ;
    TEST_MSG ("binary sum of {A,Empty,B} did not equal A+B") ;
    OK (GrB_free (&Empty)) ;
    OK (GrB_free (&C)) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    OK (GrB_free (&Expected)) ;

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Matrix_Binary_Sum_odd: exercise the odd-count carry path
//------------------------------------------------------------------------------

// Sum odd numbers of matrices (3, 5, 7, ...) so that at one or more levels an
// unpaired trailing matrix must be carried up to the next level unchanged.
// Compare against an independently accumulated expected.

void test_Matrix_Binary_Sum_odd (void)
{
    setup ( ) ;

    GrB_Index counts [ ] = { 3, 5, 7, 9 } ;
    int ncounts = sizeof (counts) / sizeof (counts [0]) ;

    for (int t = 0 ; t < ncounts ; t++)
    {
        GrB_Index nmat = counts [t] ;
        GrB_Matrix *Mats = NULL ;
        OK (LAGraph_Malloc ((void **) &Mats, nmat, sizeof (GrB_Matrix), msg)) ;

        OK (GrB_Matrix_new (&Expected, GrB_FP64, 8, 8)) ;
        for (GrB_Index k = 0 ; k < nmat ; k++)
        {
            // 3 distinct entries per matrix; (5,5) is shared by all matrices,
            // others overlap across matrices, so entries must be combined
            GrB_Index Mi [ ] = { (GrB_Index) (k % 8), 1, 5 } ;
            GrB_Index Mj [ ] = { 2, (GrB_Index) (k % 8), 5 } ;
            double    Mx [ ] = { (double) (k + 1), 1, 2 } ;
            Mats [k] = NULL ;
            OK (GrB_Matrix_new (&Mats [k], GrB_FP64, 8, 8)) ;
            OK (GrB_Matrix_build_FP64 (Mats [k], Mi, Mj, Mx, 3, NULL)) ;
            OK (GrB_eWiseAdd (Expected, NULL, NULL, GrB_PLUS_FP64, Expected,
                Mats [k], NULL)) ;
        }

        OK (LAGraph_Matrix_Binary_Sum (&C, Mats, nmat, GrB_PLUS_FP64, msg)) ;

        bool ok ;
        OK (LAGraph_Matrix_IsEqual (&ok, C, Expected, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("binary sum of %d matrices did not match expected",
            (int) nmat) ;

        for (GrB_Index k = 0 ; k < nmat ; k++)
        {
            OK (GrB_free (&Mats [k])) ;
        }
        OK (LAGraph_Free ((void **) &Mats, msg)) ;
        OK (GrB_free (&C)) ;
        OK (GrB_free (&Expected)) ;
    }

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Matrix_Binary_Sum_types: exercise every built-in type branch
//------------------------------------------------------------------------------

void test_Matrix_Binary_Sum_types (void)
{
    setup ( ) ;

    // (type, GrB_PLUS operator) pairs, one per built-in type branch
    GrB_Type   types [ ] = { GrB_BOOL, GrB_INT8, GrB_INT16, GrB_INT32,
        GrB_INT64, GrB_UINT8, GrB_UINT16, GrB_UINT32, GrB_UINT64,
        GrB_FP32, GrB_FP64 } ;
    GrB_BinaryOp ops [ ] = { GrB_LOR, GrB_PLUS_INT8, GrB_PLUS_INT16,
        GrB_PLUS_INT32, GrB_PLUS_INT64, GrB_PLUS_UINT8, GrB_PLUS_UINT16,
        GrB_PLUS_UINT32, GrB_PLUS_UINT64, GrB_PLUS_FP32, GrB_PLUS_FP64 } ;
    int ntypes = sizeof (types) / sizeof (types [0]) ;

    GrB_Index Ai [ ] = { 0, 1, 2 } ;
    GrB_Index Aj [ ] = { 0, 1, 0 } ;
    GrB_Index Bi [ ] = { 0, 1, 0 } ;
    GrB_Index Bj [ ] = { 0, 1, 2 } ;

    for (int t = 0 ; t < ntypes ; t++)
    {
        OK (GrB_Matrix_new (&A, types [t], 3, 3)) ;
        OK (GrB_Matrix_new (&B, types [t], 3, 3)) ;
        // build with INT32 values; GraphBLAS typecasts to the matrix type
        int32_t Ax [ ] = { 1, 1, 1 } ;
        int32_t Bx [ ] = { 1, 1, 1 } ;
        OK (GrB_Matrix_build_INT32 (A, Ai, Aj, Ax, 3, NULL)) ;
        OK (GrB_Matrix_build_INT32 (B, Bi, Bj, Bx, 3, NULL)) ;

        OK (GrB_Matrix_new (&Expected, types [t], 3, 3)) ;
        OK (GrB_eWiseAdd (Expected, NULL, NULL, ops [t], A, B, NULL)) ;

        GrB_Matrix Mats [2] = { A, B } ;
        OK (LAGraph_Matrix_Binary_Sum (&C, Mats, 2, ops [t], msg)) ;

        bool ok ;
        OK (LAGraph_Matrix_IsEqual (&ok, C, Expected, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("type index %d failed", t) ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
        OK (GrB_free (&C)) ;
        OK (GrB_free (&Expected)) ;
    }

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Matrix_Binary_Sum_brutal
//------------------------------------------------------------------------------

#if LG_BRUTAL_TESTS
void test_Matrix_Binary_Sum_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    // five matrices: an odd count, so the carry path is exercised under brutal
    // malloc-failure testing as well, validating the Pool / W / Wnext cleanup
    #define NBRUTAL 5
    GrB_Matrix Mats [NBRUTAL] ;
    OK (GrB_Matrix_new (&Expected, GrB_FP64, 6, 6)) ;
    for (int k = 0 ; k < NBRUTAL ; k++)
    {
        GrB_Index Mi [ ] = { (GrB_Index) (k % 6), 1, 4 } ;
        GrB_Index Mj [ ] = { 2, (GrB_Index) (k % 6), 4 } ;
        double    Mx [ ] = { (double) (k + 1), 1, 2 } ;
        Mats [k] = NULL ;
        OK (GrB_Matrix_new (&Mats [k], GrB_FP64, 6, 6)) ;
        OK (GrB_Matrix_build_FP64 (Mats [k], Mi, Mj, Mx, 3, NULL)) ;
        OK (GrB_eWiseAdd (Expected, NULL, NULL, GrB_PLUS_FP64, Expected,
            Mats [k], NULL)) ;
    }

    LG_BRUTAL (LAGraph_Matrix_Binary_Sum (&C, Mats, NBRUTAL, GrB_PLUS_FP64,
        msg)) ;

    bool ok ;
    OK (LAGraph_Matrix_IsEqual (&ok, C, Expected, msg)) ;
    TEST_CHECK (ok) ;

    for (int k = 0 ; k < NBRUTAL ; k++)
    {
        OK (GrB_free (&Mats [k])) ;
    }
    OK (GrB_free (&C)) ;
    OK (GrB_free (&Expected)) ;
    #undef NBRUTAL

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//------------------------------------------------------------------------------
// test_Matrix_Binary_Sum_parallel: exercise the parallel pair-sum path
//------------------------------------------------------------------------------

// Sum many matrices with multiple outer threads and confirm the result matches
// an independently accumulated expected.  This stresses the concurrent
// pair-sums of the binary reduction under real outer parallelism, including an
// odd count to also drive the carry path.

void test_Matrix_Binary_Sum_parallel (void)
{
    setup ( ) ;

    // request 4 outer threads (saving and restoring the prior settings)
    int save_outer, save_inner ;
    OK (LAGraph_GetNumThreads (&save_outer, &save_inner, msg)) ;
    OK (LAGraph_SetNumThreads (4, save_inner, msg)) ;

    GrB_Index counts [ ] = { 16, 15 } ;     // even and odd
    int ncounts = sizeof (counts) / sizeof (counts [0]) ;

    for (int t = 0 ; t < ncounts ; t++)
    {
        GrB_Index nmat = counts [t] ;
        GrB_Matrix *Mats = NULL ;
        OK (LAGraph_Malloc ((void **) &Mats, nmat, sizeof (GrB_Matrix), msg)) ;

        OK (GrB_Matrix_new (&Expected, GrB_FP64, 10, 10)) ;
        for (GrB_Index k = 0 ; k < nmat ; k++)
        {
            // each matrix has 3 distinct (i,j) entries; the (7,7) entry is
            // shared by every matrix and others overlap across matrices, so
            // duplicates must be combined when the matrices are summed
            GrB_Index Mi [ ] = { (GrB_Index) (k % 10), 2, 7 } ;
            GrB_Index Mj [ ] = { 3, (GrB_Index) (k % 10), 7 } ;
            double    Mx [ ] = { (double) (k + 1), 1, 2 } ;
            Mats [k] = NULL ;
            OK (GrB_Matrix_new (&Mats [k], GrB_FP64, 10, 10)) ;
            OK (GrB_Matrix_build_FP64 (Mats [k], Mi, Mj, Mx, 3, NULL)) ;
            // accumulate into Expected independently
            OK (GrB_eWiseAdd (Expected, NULL, NULL, GrB_PLUS_FP64, Expected,
                Mats [k], NULL)) ;
        }

        OK (LAGraph_Matrix_Binary_Sum (&C, Mats, nmat, GrB_PLUS_FP64, msg)) ;

        bool ok ;
        OK (LAGraph_Matrix_IsEqual (&ok, C, Expected, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("parallel binary sum of %d matrices did not match expected",
            (int) nmat) ;

        for (GrB_Index k = 0 ; k < nmat ; k++)
        {
            OK (GrB_free (&Mats [k])) ;
        }
        OK (LAGraph_Free ((void **) &Mats, msg)) ;
        OK (GrB_free (&C)) ;
        OK (GrB_free (&Expected)) ;
    }

    // restore the original thread settings
    OK (LAGraph_SetNumThreads (save_outer, save_inner, msg)) ;

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Matrix_Binary_Sum_failures: test error handling
//------------------------------------------------------------------------------

void test_Matrix_Binary_Sum_failures (void)
{
    setup ( ) ;

    make_AB ( ) ;
    GrB_Matrix Mats [2] = { A, B } ;
    int result ;

    // NULL output
    result = LAGraph_Matrix_Binary_Sum (NULL, Mats, 2, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // NULL Matrices array
    C = NULL ;
    result = LAGraph_Matrix_Binary_Sum (&C, NULL, 2, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;
    TEST_CHECK (C == NULL) ;

    // nmatrices == 0
    result = LAGraph_Matrix_Binary_Sum (&C, Mats, 0, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_INVALID_VALUE) ;

    // NULL dup operator
    result = LAGraph_Matrix_Binary_Sum (&C, Mats, 2, NULL, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // NULL entry in the array
    GrB_Matrix WithNull [2] = { A, NULL } ;
    result = LAGraph_Matrix_Binary_Sum (&C, WithNull, 2, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // dimension mismatch
    GrB_Matrix D = NULL ;
    OK (GrB_Matrix_new (&D, GrB_FP64, 5, 5)) ;
    GrB_Matrix BadDim [2] = { A, D } ;
    result = LAGraph_Matrix_Binary_Sum (&C, BadDim, 2, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_DIMENSION_MISMATCH) ;
    OK (GrB_free (&D)) ;

    // type mismatch
    GrB_Matrix E = NULL ;
    OK (GrB_Matrix_new (&E, GrB_INT32, 4, 4)) ;
    GrB_Matrix BadType [2] = { A, E } ;
    result = LAGraph_Matrix_Binary_Sum (&C, BadType, 2, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_DOMAIN_MISMATCH) ;
    OK (GrB_free (&E)) ;

    // user-defined type not supported
    typedef struct { double x, y ; } udt_t ;
    GrB_Type UDT = NULL ;
    OK (GrB_Type_new (&UDT, sizeof (udt_t))) ;
    GrB_Matrix U = NULL ;
    OK (GrB_Matrix_new (&U, UDT, 4, 4)) ;
    GrB_Matrix Udt [1] = { U } ;
    result = LAGraph_Matrix_Binary_Sum (&C, Udt, 1, GrB_PLUS_FP64, msg) ;
    TEST_CHECK (result == GrB_NOT_IMPLEMENTED) ;
    OK (GrB_free (&U)) ;
    OK (GrB_free (&UDT)) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//------------------------------------------------------------------------------

TEST_LIST =
{
    { "Matrix_Binary_Sum", test_Matrix_Binary_Sum },
    { "Matrix_Binary_Sum_odd", test_Matrix_Binary_Sum_odd },
    { "Matrix_Binary_Sum_types", test_Matrix_Binary_Sum_types },
    { "Matrix_Binary_Sum_parallel", test_Matrix_Binary_Sum_parallel },
    { "Matrix_Binary_Sum_failures", test_Matrix_Binary_Sum_failures },
    #if LG_BRUTAL_TESTS
    { "Matrix_Binary_Sum_brutal", test_Matrix_Binary_Sum_brutal },
    #endif
    { NULL, NULL }
} ;
