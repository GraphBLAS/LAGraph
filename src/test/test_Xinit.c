//------------------------------------------------------------------------------
// LAGraph/src/test/test_Xinit.c:  test LAGraph_Xinit and LAGraph_Global
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
// test_Xinit:  test LAGraph_Xinit
//------------------------------------------------------------------------------

void test_Xinit (void)
{

    printf ("\nTesting LAGraph_Xinit:\n") ;

    TEST_CHECK (LAGraph_Xinit (NULL, NULL, NULL, NULL, true, msg) == -1) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_Xinit (malloc, NULL, NULL, NULL, true, msg) == -1) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_Xinit (NULL, NULL, NULL, free, true, msg) == -1) ;
    printf ("msg: %s\n", msg) ;

    OK (LAGraph_Xinit (malloc, NULL, NULL, free, true, msg)) ;
    printf ("msg: [%s]\n", msg) ;

    // LAGraph_Xinit cannot be called twice
    int status = LAGraph_Xinit (malloc, calloc, realloc, free, true, msg) ;
    TEST_CHECK (status == -GrB_INVALID_VALUE) ;

    // TODO: this error message is not informative
    printf ("msg: %s\n", msg) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Xinit", test_Xinit },
    { NULL, NULL }
} ;

