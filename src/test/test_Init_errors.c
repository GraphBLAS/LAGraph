//------------------------------------------------------------------------------
// LAGraph/src/test/test_Init_errors.c:  test LAGraph_Init and LAGraph_Finalize
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LAGraph_test.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test_Init_errors:  test LAGraph_Init
//------------------------------------------------------------------------------

void test_Init_errors (void)
{

    int status = LAGraph_Init (msg) ;
    OK (status) ;
    int ver [3] ;

    // LAGraph_Init cannot be called twice
    status = LAGraph_Init (msg) ;
    printf ("\nstatus: %d msg: %s\n", status, msg) ;
    TEST_CHECK (status != GrB_SUCCESS) ;

    OK (LAGraph_Finalize (msg)) ;

    // calling LAGraph_Finalize twice leads to undefined behavior;
    // for SuiteSparse, it returns GrB_SUCCESS
    status = LAGraph_Finalize (msg) ;
    printf ("status %d\n", status) ;
    #if LAGRAPH_SUITESPARSE
    TEST_CHECK (status == GrB_SUCCESS) ;
    #endif
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Init_errors", test_Init_errors },
    // no brutal test: see test_Xinit
    { NULL, NULL }
} ;

