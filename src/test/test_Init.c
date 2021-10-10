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

    #if LG_SUITESPARSE
    const char *name, *date ;
    int ver [3] ;
    OK (GxB_get (GxB_LIBRARY_NAME, &name)) ;
    OK (GxB_get (GxB_LIBRARY_DATE, &date)) ;
    OK (GxB_get (GxB_LIBRARY_VERSION, ver)) ;
    printf ("\n%s %d.%d.%d (%s)\n", name, ver [0], ver [1], ver [2], date) ;
    #else
    printf ("\nVanilla GraphBLAS: no GxB* extensions\n") ;
    #endif

    // LAGraph_Init cannot be called twice
    TEST_CHECK (LAGraph_Init (msg) != GrB_SUCCESS) ;

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

