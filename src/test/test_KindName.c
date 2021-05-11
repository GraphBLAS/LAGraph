//------------------------------------------------------------------------------
// LAGraph/src/test/test_KindName.c:  test LAGraph_KindName
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
char *name = NULL ;

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

    OK (LAGraph_KindName (&name, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
    OK (strcmp (name, "undirected")) ;

    OK (LAGraph_KindName (&name, LAGRAPH_ADJACENCY_DIRECTED, msg)) ;
    OK (strcmp (name, "directed")) ;

    OK (LAGraph_KindName (&name, LAGRAPH_UNKNOWN, msg)) ;
    OK (strcmp (name, "unknown")) ;

    TEST_CHECK (LAGraph_KindName (&name, 42, msg) == -1) ;
    printf ("\nmsg: %s\n", msg) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "KindName", test_KindName },
    { NULL, NULL }
} ;

