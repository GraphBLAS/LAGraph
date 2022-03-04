//------------------------------------------------------------------------------
// LAGraph/src/test/test_Init.c:  test LAGraph_Init and LAGraph_Finalize
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
// test_Init:  test LAGraph_Init
//------------------------------------------------------------------------------

void test_Init (void)
{

    int status = LAGraph_Init (msg) ;
    OK (status) ;
    int ver [3] ;

    #if LAGRAPH_SUITESPARSE
    const char *name, *date ;
    OK (GxB_get (GxB_LIBRARY_NAME, &name)) ;
    OK (GxB_get (GxB_LIBRARY_DATE, &date)) ;
    OK (GxB_get (GxB_LIBRARY_VERSION, ver)) ;
    printf ("\nlibrary: %s %d.%d.%d (%s)\n", name, ver [0], ver [1], ver [2],
        date) ;
    printf (  "include: %s %d.%d.%d (%s)\n", GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR, GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB, GxB_IMPLEMENTATION_DATE) ;
    // make sure the SuiteSparse:GraphBLAS version and date match
    TEST_CHECK (ver [0] == GxB_IMPLEMENTATION_MAJOR) ;
    TEST_CHECK (ver [1] == GxB_IMPLEMENTATION_MINOR) ;
    TEST_CHECK (ver [2] == GxB_IMPLEMENTATION_SUB) ;
    OK (strcmp (date, GxB_IMPLEMENTATION_DATE)) ;
    #else
    printf ("\nVanilla GraphBLAS: no GxB* extensions\n") ;
    #endif

    // check the LAGraph version using both LAGraph.h and LAGraph_Version
    printf ("LAGraph version %d.%d.%d (%s) from LAGraph.h\n",
        LAGRAPH_VERSION_MAJOR, LAGRAPH_VERSION_MINOR, LAGRAPH_VERSION_UPDATE,
        LAGRAPH_DATE) ;

    char version_date [LAGRAPH_MSG_LEN] ;
    status = LAGraph_Version (ver, version_date, msg) ;
    OK (status) ;

    printf ("LAGraph version %d.%d.%d (%s) from LAGraph_Version\n",
        ver [0], ver [1], ver [2], version_date) ;

    // make sure the LAGraph version and date match
    TEST_CHECK (ver [0] == LAGRAPH_VERSION_MAJOR) ;
    TEST_CHECK (ver [1] == LAGRAPH_VERSION_MINOR) ;
    TEST_CHECK (ver [2] == LAGRAPH_VERSION_UPDATE) ;
    OK (strcmp (version_date, LAGRAPH_DATE)) ;

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
    { "Init", test_Init },
    // no brutal test: see test_Xinit
    { NULL, NULL }
} ;

