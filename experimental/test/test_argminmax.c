//------------------------------------------------------------------------------
// experimental/test/test_argminmax:  test LAGraph_argminmax
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraph.h>
#include <LG_internal.h>

char msg [LAGRAPH_MSG_LEN] ;

#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    { "structure.mtx" },
    { "karate.mtx" },
    { "west0067.mtx" },
    { "bcsstk13.mtx" },
    { "" },
} ;

void test_argminmax (void)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, C = NULL ;
    GrB_Vector x = NULL, p = NULL, x2 = NULL, p2 = NULL ;
    GrB_Index nrows, ncols ;

    for (int k = 0 ; ; k++)
    {
        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break ;
        // printf ("\n %s: ==================================\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;
        fclose (f) ;
        OK (GrB_Matrix_nrows (&nrows, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;

        for (int ktype = 0 ; ktype < 11 ; ktype++)
        {
            GrB_Type type ;
            switch (ktype)
            {
                case  0: type = GrB_BOOL    ; break ;
                case  1: type = GrB_INT8    ; break ;
                case  2: type = GrB_INT16   ; break ;
                case  3: type = GrB_INT32   ; break ;
                case  4: type = GrB_INT64   ; break ;
                case  5: type = GrB_UINT8   ; break ;
                case  6: type = GrB_UINT16  ; break ;
                case  7: type = GrB_UINT32  ; break ;
                case  8: type = GrB_UINT64  ; break ;
                case  9: type = GrB_FP32    ; break ;
                default:
                case 10: type = GrB_FP64    ; break ;
            }

            // typecast A into a different type
            OK (GrB_Matrix_new (&C, type, nrows, ncols)) ;
            OK (GrB_assign (C, NULL, NULL, A,
                GrB_ALL, nrows, GrB_ALL, ncols, NULL)) ;

            // printf ("\nA:\n") ;
            // OK (LAGraph_Matrix_Print (A, 2, stdout, msg)) ;

            // printf ("\nC:\n") ;
            // OK (LAGraph_Matrix_Print (C, 2, stdout, msg)) ;

            for (int is_min = 0 ; is_min <= 1 ; is_min++)
            {
                for (int dim = 0 ; dim <= 2 ; dim++)
                {
                    // printf ("\nis_min: %d dim: %d\n", is_min, dim) ;
                    // test the algorithm
                    OK (LAGraph_argminmax (&x, &p, C, dim, is_min, msg)) ;
                    // printf ("\nx:\n") ;
                    // OK (LAGraph_Vector_Print (x, 2, stdout, msg)) ;
                    // printf ("\np:\n") ;
                    // OK (LAGraph_Vector_Print (p, 2, stdout, msg)) ;
                    // check the result
                    OK (LG_check_argminmax (&x2, &p2, C, dim, is_min, msg)) ;
                    // printf ("\nx2:\n") ;
                    // OK (LAGraph_Vector_Print (x2, 2, stdout, msg)) ;
                    // printf ("\np2:\n") ;
                    // OK (LAGraph_Vector_Print (p2, 2, stdout, msg)) ;
                    bool isequal = false ;
                    // x and x2 must be equal, for all cases
                    OK (LAGraph_Vector_IsEqual (&isequal, x, x2, msg)) ;
                    TEST_CHECK (isequal) ;
                    uint64_t npvals = 0 ;
                    OK (GrB_Vector_nvals (&npvals, p)) ;
                    if (dim > 0 || npvals == 0)
                    {
                        // For dim=1 or dim=2 (row-wise or col-wise), p and p2
                        // must always match
                        OK (LAGraph_Vector_IsEqual (&isequal, p, p2, msg)) ;
                        TEST_CHECK (isequal) ;
                    }
                    else
                    {
                        // For dim=0, the result is a single scalar, with
                        // C(p2[0],p2[1]) being argmin/argmax of C.  The two
                        // methods may find different places where the min/max
                        // entry appears in C, if there are ties, so p and p2
                        // can differ.  Just make sure C(p2[0],p2[1]) is equal
                        // to x [0].
                        uint64_t i, j ;
                        TEST_CHECK (npvals == 2) ;
                        GrB_Info info = GrB_Vector_extractElement (&i, p2, 0) ;
                        TEST_CHECK (info >= GrB_SUCCESS) ;
                        info = GrB_Vector_extractElement (&j, p2, 1) ;
                        TEST_CHECK (info >= GrB_SUCCESS) ;
                        double x_1 = 0, x_2 = 1 ;
                        info = GrB_Matrix_extractElement (&x_1, C, i, j) ;
                        TEST_CHECK (info >= GrB_SUCCESS) ;
                        info = GrB_Vector_extractElement (&x_2, x, 0) ;
                        TEST_CHECK (info >= GrB_SUCCESS) ;
                        // printf ("x_1 %g x_2 %g\n", x_1, x_2) ;
                        TEST_CHECK (x_1 == x_2) ;
                    }
                    OK (GrB_free (&x)) ;
                    OK (GrB_free (&p)) ;
                    OK (GrB_free (&x2)) ;
                    OK (GrB_free (&p2)) ;
                }
            }
            OK (GrB_free (&C)) ;
        }
        OK (GrB_free (&A)) ;
    }

    //--------------------------------------------------------------------------
    // finalize LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Finalize (msg) ;
#endif
}

//----------------------------------------------------------------------------
// test_argminmax_errors
//----------------------------------------------------------------------------

void test_argminmax_errors (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Vector x = NULL, p = NULL ;
    OK (LAGraph_Random_Matrix (&A, GrB_FP64, 5, 5, 0.5, 1, msg)) ;
    GrB_Info info = LG_check_argminmax (&x, &p, A, 3, true, msg) ;
    TEST_CHECK (info == GrB_INVALID_VALUE) ;
    GrB_free (&A) ;
    LAGraph_Finalize (msg) ;
#endif
}

//----------------------------------------------------------------------------
// the main program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"argminmax", test_argminmax},
    {"argminmax_errors", test_argminmax_errors},
    {NULL, NULL}
} ;

