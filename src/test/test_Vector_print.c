//-----------------------------------------------------------------------------
// LAGraph/src/test/test_Vector_print.c: test LAGraph_Vector_print
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include "LAGraph_test.h"

//-----------------------------------------------------------------------------
// test_print
//-----------------------------------------------------------------------------

void test_print (void)
{
    char msg [LAGRAPH_MSG_LEN] ;
    OK (LAGraph_Init (msg)) ;

    GrB_Vector v = NULL ;
    GrB_Index n = 40 ;

    OK (GrB_Vector_new (&v, GrB_BOOL, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 1, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_INT8, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, -8, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_INT16, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, -16, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_INT32, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, -32, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_INT64, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, -64, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_UINT8, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 8, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_UINT16, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 16, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_UINT32, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 32, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_UINT64, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 64, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_FP32, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 3.14159, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Vector_new (&v, GrB_FP64, n)) ;
    OK (GrB_assign (v, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    OK (GrB_Vector_setElement (v, 99.999, 0)) ;
    OK (LAGraph_Vector_print (v, 2, stdout, msg)) ;

    // attempt to print to a NULL file, which should fail
    int result = LAGraph_Vector_print (v, 2, NULL, msg) ;
    TEST_CHECK (result == -1001) ;
    OK (GrB_Vector_free (&v)) ;

    // attempt to print a vector with a user-defined type, which should fail
    GrB_Type type = NULL ;
    OK (GrB_Type_new (&type, sizeof (int))) ;
    OK (GrB_Vector_new (&v, type, n)) ;
    result = LAGraph_Vector_print (v, 2, stdout, msg) ;
    TEST_CHECK (result == -1002) ;
    OK (GrB_Vector_free (&v)) ;

    OK (GrB_Type_free (&type)) ;
    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST = {
    {"test_print", test_print},
    {NULL, NULL}
};
