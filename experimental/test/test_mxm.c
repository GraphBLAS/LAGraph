//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_mxm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>

//------------------------------------------------------------------------------
// test cases
//------------------------------------------------------------------------------

char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL, B = NULL, C2 = NULL, Cin = NULL, C = NULL ;

//------------------------------------------------------------------------------
// test_mxm: test GrB_mxm with different types of inputs
//------------------------------------------------------------------------------

void test_mxm (void)
{

    GrB_Info info ;
    OK (LAGraph_Init (msg)) ;
    uint64_t n = 256 ;
    uint64_t seed = 0 ;

    for (int t = 32 ; t <= 64 ; t += 32)
    {
        // set up the type and semiring
        printf ("\ntest_mxm ======================= type : FP%d\n",  t) ;
        GrB_Type type = (t == 32) ? GrB_FP32 : GrB_FP64 ;
        GrB_Semiring semiring = (t == 32) ?
            GrB_PLUS_TIMES_SEMIRING_FP32 : GrB_PLUS_TIMES_SEMIRING_FP64 ;
        GrB_BinaryOp accum = (t == 32) ?
            GrB_PLUS_FP32 : GrB_PLUS_FP64 ;
        GrB_BinaryOp minus = (t == 32) ?
            GrB_MINUS_FP32 : GrB_MINUS_FP64 ;
        GrB_UnaryOp abs_op = (t == 32) ?
            GrB_ABS_FP32 : GrB_ABS_FP64 ;
        GrB_Monoid max_monoid = (t == 32) ?
            GrB_MAX_MONOID_FP32 : GrB_MAX_MONOID_FP64 ;
        double tol = (t == 32) ? 1e-4 : 1e-11 ;

        // create some random test matrices
        OK (LAGraph_Random_Matrix (&A, type, n, n, INFINITY, seed, msg)) ;
        seed += n*n ;
        OK (LAGraph_Random_Matrix (&B, type, n, n, INFINITY, seed, msg)) ;
        seed += n*n ;
        OK (LAGraph_Random_Matrix (&Cin, type, n, n, INFINITY, seed, msg)) ;
        seed += n*n ;

        // C = Cin + A*B
        OK (GrB_Matrix_dup (&C, Cin)) ;
        OK (GrB_mxm (C, NULL, accum, semiring, A, B, NULL)) ;
        double maxerr = 0 ;

        // test with different sparsity formats and with JIT on/off
        for (int A_sparsity = 1 ; A_sparsity <= 8 ; A_sparsity *= 2)
        {
            LG_SET_FORMAT_HINT (A, A_sparsity) ;
            for (int B_sparsity = 1 ; B_sparsity <= 8 ; B_sparsity *= 2)
            {
                LG_SET_FORMAT_HINT (B, B_sparsity) ;
                for (int C_sparsity = 1 ; C_sparsity <= 8 ; C_sparsity *= 2)
                {
                    LG_SET_FORMAT_HINT (Cin, C_sparsity) ;
                    for (int jit = 0 ; jit <= 1 ; jit++)
                    {
                        OK (LG_SET_JIT (jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;

                        // C2 = Cin + A*B
                        OK (GrB_Matrix_dup (&C2, Cin)) ;
                        OK (GrB_mxm (C2, NULL, accum, semiring, A, B, NULL)) ;

                        // C2 = abs (C - C2)
                        OK (GrB_eWiseAdd (C2, NULL, NULL, minus, C2, C, NULL)) ;
                        OK (GrB_apply (C2, NULL, NULL, abs_op, C2, NULL)) ;

                        // err = max (C2)
                        double err = 0 ;
                        OK (GrB_reduce (&err, NULL, max_monoid, C2, NULL)) ;
                        TEST_CHECK (err < tol) ;
                        maxerr = LAGRAPH_MAX (maxerr, err) ;
                        GrB_free (&C2) ;
                    }
                }
            }
        }

        printf ("max err: %g\n", maxerr) ;
        GrB_free (&A) ;
        GrB_free (&B) ;
        GrB_free (&C) ;
        GrB_free (&Cin) ;
    }

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "mxm", test_mxm },
    { NULL, NULL }
} ;

