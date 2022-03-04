//-----------------------------------------------------------------------------
// LAGraph/src/test/test_Multiply_size_t.c: test LG_Multiply_size_t
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//-----------------------------------------------------------------------------

#include "LAGraph_test.h"
#include "LG_internal.h"

//-----------------------------------------------------------------------------
// test_multiply
//-----------------------------------------------------------------------------

void test_multiply (void)
{
    char msg [LAGRAPH_MSG_LEN] ;
    OK (LAGraph_Init (msg)) ;

    size_t c = 99 ;
    TEST_CHECK (LG_Multiply_size_t (&c, (size_t) 0, (size_t) 42)) ;
    TEST_CHECK (c == 0) ;

    TEST_CHECK (LG_Multiply_size_t (&c, (size_t) 77, (size_t) 42)) ;
    TEST_CHECK (c == 77 * 42) ;

    TEST_CHECK (!LG_Multiply_size_t (&c, (size_t) SIZE_MAX, (size_t) 42)) ;

    size_t a = SIZE_MAX / 2 ;
    TEST_CHECK (!LG_Multiply_size_t (&c, a, a)) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST = {
    {"test_multiply", test_multiply},
    // no brutal test needed
    {NULL, NULL}
};
