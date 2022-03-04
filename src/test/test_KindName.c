//------------------------------------------------------------------------------
// LAGraph/src/test/test_KindName.c:  test LAGraph_KindName
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
char name [LAGRAPH_MAX_NAME_LEN] ;

//------------------------------------------------------------------------------
// setup: start a test
//------------------------------------------------------------------------------

void setup (void)
{
    OK (LAGraph_Init (msg)) ;
}

//------------------------------------------------------------------------------
// teardown: finalize a test
//------------------------------------------------------------------------------

void teardown (void)
{
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_KindName:  test LAGraph_KindName
//------------------------------------------------------------------------------

void test_KindName (void)
{
    setup ( ) ;

    OK (LAGraph_KindName (name, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    OK (strcmp (name, "undirected")) ;

    OK (LAGraph_KindName (name, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    OK (strcmp (name, "directed")) ;

    OK (LAGraph_KindName (name, LAGRAPH_UNKNOWN, msg)) ;
    OK (strcmp (name, "unknown")) ;

    TEST_CHECK (LAGraph_KindName (name, 42, msg) == GrB_INVALID_VALUE) ;
    printf ("\nmsg: %s\n", msg) ;

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_KindName_brutal
//------------------------------------------------------------------------------

// LAGraph_KindName currently doesn't do any mallocs so this test is not
// strictly necessary, but it's simple to include here.  It serves as a very
// simple use-case of the brutal testing mechanism.

#if LAGRAPH_SUITESPARSE
void test_KindName_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;
    LG_BRUTAL (LAGraph_KindName (name, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    OK (strcmp (name, "undirected")) ;
    OK (LG_brutal_teardown (msg)) ;
}
#endif

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "KindName", test_KindName },
    { "KindName_brutal", test_KindName_brutal },
    { NULL, NULL }
} ;

