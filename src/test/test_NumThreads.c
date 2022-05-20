//------------------------------------------------------------------------------
// LAGraph/src/test/test_NumThreads.c:  test LAGraph_(Get,Set)NumThreads
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

int nthreads_hi = 0, nthreads_lo = 0 ;
char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test_NumThreads:  test LAGraph_GetNumThreads and LAGraph_SetNumThreads
//------------------------------------------------------------------------------

void test_NumThreads (void)
{

    OK (LAGraph_Init (msg)) ;

    nthreads_hi = 0 ;
    nthreads_lo = 0 ;
    OK (LAGraph_GetNumThreads (&nthreads_hi, &nthreads_lo, msg)) ;
    TEST_CHECK (nthreads_hi > 0) ;
    TEST_CHECK (nthreads_lo > 0) ;

    nthreads_hi = 0 ;
    nthreads_lo = 0 ;
    OK (LAGraph_GetNumThreads (&nthreads_hi, &nthreads_lo, NULL)) ;
    TEST_CHECK (nthreads_hi > 0) ;
    TEST_CHECK (nthreads_lo > 0) ;

    OK (LAGraph_SetNumThreads (2, 4, msg)) ;
    nthreads_hi = 0 ;
    nthreads_lo = 0 ;
    OK (LAGraph_GetNumThreads (&nthreads_hi, &nthreads_lo, msg)) ;
    TEST_CHECK (nthreads_hi > 0) ;
    TEST_CHECK (nthreads_lo > 0) ;

    OK (LAGraph_SetNumThreads (2, 4, NULL)) ;
    nthreads_hi = 0 ;
    nthreads_lo = 0 ;
    OK (LAGraph_GetNumThreads (&nthreads_hi, &nthreads_lo, NULL)) ;
    TEST_CHECK (nthreads_hi > 0) ;
    TEST_CHECK (nthreads_lo > 0) ;

    TEST_CHECK (LAGraph_GetNumThreads (NULL, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("\nmsg: %s\n", msg) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "NumThreads", test_NumThreads },
    // no brutal test needed
    { NULL, NULL }
} ;

