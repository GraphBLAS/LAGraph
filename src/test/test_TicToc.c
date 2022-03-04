//------------------------------------------------------------------------------
// LAGraph/src/test/test_TicToc.c:  test LAGraph_Tic and LAGraph_Toc
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

double t, tic [2] ;
char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test_TicToc:  test LAGraph_TicToc
//------------------------------------------------------------------------------

void test_TicToc (void)
{

    OK (LAGraph_Init (msg)) ;

    // start the timer
    OK (LAGraph_Tic (tic, msg)) ;

    // do some useless work
    double x = msg [0] ;
    for (int64_t k = 0 ; k < 10000 ; k++)
    {
        for (int64_t i = 0 ; i < 10000 ; i++)
        {
            x = x + 1 ;
            if (x > 100) x = x/2 ;
        }
    }

    // stop the timer
    OK (LAGraph_Toc (&t, tic, msg)) ;

    // print the result so the compiler doesn't remove the loops above
    printf ("\nresult: %g, time: %g sec\n", x, t) ;

    // check error conditions
    TEST_CHECK (LAGraph_Toc (NULL, NULL, msg) == GrB_NULL_POINTER) ;
    TEST_CHECK (LAGraph_Toc (&t,   NULL, msg) == GrB_NULL_POINTER) ;
    TEST_CHECK (LAGraph_Toc (NULL, tic,  msg) == GrB_NULL_POINTER) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "TicToc", test_TicToc },
    // no brutal test needed
    { NULL, NULL }
} ;

