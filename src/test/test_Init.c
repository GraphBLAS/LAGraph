//------------------------------------------------------------------------------
// LAGraph/src/test/test_Init.c:  test LAGraph_Init and LAGraph_Finalize
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LAGraph_test.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test_Init:  test LAGraph_Init
//------------------------------------------------------------------------------

void test_Init (void)
{

    OK (LAGraph_Init (msg)) ;

    // LAGraph_Init cannot be called twice
    TEST_CHECK (LAGraph_Init (msg) == -GrB_INVALID_VALUE) ;

    // TODO: this error message is not informative
    printf ("\nmsg: %s\n", msg) ;

    OK (LAGraph_Finalize (msg)) ;

    // calling LAGraph_Finalize twice leads to undefined behavior;
    // for SuiteSparse, it returns GrB_SUCCESS
    int status = LAGraph_Finalize (msg) ;
    printf ("status %d\n", status) ;
    #if LG_SUITESPARSE
    TEST_CHECK (status == GrB_SUCCESS) ;
    #endif
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Init", test_Init },
    { NULL, NULL }
} ;

