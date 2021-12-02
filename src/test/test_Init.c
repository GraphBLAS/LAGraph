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
    printf ("\nlibrary: %s %d.%d.%d (%s)\n", name, ver [0], ver [1], ver [2],
        date) ;
    printf (  "include: %s %d.%d.%d (%s)\n", GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR, GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB, GxB_IMPLEMENTATION_DATE) ;
    #else
    printf ("\nVanilla GraphBLAS: no GxB* extensions\n") ;
    #endif

    // LAGraph_Init cannot be called twice
    int status = LAGraph_Init (msg) ;
    printf ("\nstatus: %d msg: %s\n", status, msg) ;
    TEST_CHECK (status != GrB_SUCCESS) ;

    OK (LAGraph_Finalize (msg)) ;

    // calling LAGraph_Finalize twice leads to undefined behavior;
    // for SuiteSparse, it returns GrB_SUCCESS
    status = LAGraph_Finalize (msg) ;
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
    // no brutal test: see test_Xinit
    { NULL, NULL }
} ;

