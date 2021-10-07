//------------------------------------------------------------------------------
// LAGraph/src/test/test_Random.cpp: test cases for random vector generator
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
GrB_Vector Seed = NULL ;
GrB_Vector X = NULL ;
GrB_Vector I = NULL ;

//------------------------------------------------------------------------------
// test_Random
//------------------------------------------------------------------------------

void test_Random (void)
{
    LAGraph_Init (msg) ;
    OK (LAGraph_Random_Init (msg)) ;

    int64_t seed = 42 ;

    for (int trial = 0 ; trial <= 1 ; trial++)
    {
        printf ("\n================================ seed: %ld\n", seed) ;

        // generate a seed vector (all entries present)
        printf ("\nDense random seed:\n") ;
        GrB_Index n = 8 ;
        OK (GrB_Vector_new (&Seed, GrB_INT64, n)) ;
        OK (GrB_Vector_assign_INT64 (Seed, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
        OK (LAGraph_Random_Seed (Seed, seed, msg)) ;
        OK (GxB_Vector_fprint (Seed, "Seed (dense)", GxB_COMPLETE, NULL)) ;

        // generate a random int64 vector
        printf ("\nDense random int64:\n") ;
        OK (GrB_Vector_new (&I, GrB_INT64, n)) ;
        OK (LAGraph_Random_INT64 (I, Seed, msg)) ;
        OK (GxB_Vector_fprint (I, "I (dense)", GxB_COMPLETE, NULL)) ;

        // generate a random double vector
        printf ("\nDense random double:\n") ;
        OK (GrB_Vector_new (&X, GrB_FP64, n)) ;
        OK (LAGraph_Random_FP64 (X, Seed, msg)) ;
        OK (GxB_Vector_fprint (X, "X (dense)", GxB_COMPLETE, NULL)) ;

        // free all workspace
        OK (GrB_Vector_free (&Seed)) ;
        OK (GrB_Vector_free (&X)) ;
        OK (GrB_Vector_free (&I)) ;

        // generate a sparse seed vector
        printf ("\nSparse random seed:\n") ;
        OK (GrB_Vector_new (&Seed, GrB_INT64, n)) ;
        for (int i = 0 ; i < n ; i += 2)
        {
            OK (GrB_Vector_setElement_INT64 (Seed, 0, i)) ;
        }
        OK (LAGraph_Random_Seed (Seed, seed, msg)) ;
        OK (GxB_Vector_fprint (Seed, "Seed (sparse)", GxB_COMPLETE, NULL)) ;

        // generate a random int64 vector
        printf ("\nSparse random int64:\n") ;
        OK (GrB_Vector_new (&I, GrB_INT64, n)) ;
        OK (LAGraph_Random_INT64 (I, Seed, msg)) ;
        OK (GxB_Vector_fprint (I, "I (sparse)", GxB_COMPLETE, NULL)) ;

        // generate a random double vector
        printf ("\nSparse random double: n %lu\n", n) ;
        OK (GrB_Vector_new (&X, GrB_FP64, n)) ;
        OK (LAGraph_Random_FP64 (X, Seed, msg)) ;
        OK (GxB_Vector_fprint (X, "X (sparse)", GxB_COMPLETE, NULL)) ;

        // free all workspace
        OK (GrB_Vector_free (&Seed)) ;
        OK (GrB_Vector_free (&X)) ;
        OK (GrB_Vector_free (&I)) ;

        seed = 987 ;
    }

    OK (LAGraph_Random_Finalize (msg)) ;
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// Test list
//------------------------------------------------------------------------------

TEST_LIST = {
    {"Random", test_Random},
    {NULL, NULL}
};

