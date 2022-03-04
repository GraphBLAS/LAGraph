//-----------------------------------------------------------------------------
// LAGraph/src/test/test_acutest.c: simple demo of how to use acutest
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Scott McMillan, SEI Carnegie Mellon University

//-----------------------------------------------------------------------------

#include "LAGraph_test.h"

//****************************************************************************
void test_dummy(void)
{
    TEST_CHECK(42 == 42);
    // this test message will not appear in the log, because the test passes:
    TEST_MSG ("Testing equality %d", 42) ;
}

//****************************************************************************
#if 0
void test_dummy_fails(void)
{
    TEST_CHECK(42 == 0);
    // this test message will appear in the log, because the test fails:
    TEST_MSG ("this test is supposed to fail, because %d != %d\n", 42, 0) ;
}
#endif

//****************************************************************************
//****************************************************************************

// uncomment the line below to create an intentional test failure.

TEST_LIST = {
    {"dummy", test_dummy},
//  {"dummy_fails", test_dummy_fails},
    // no brutal test needed
    {NULL, NULL}
};
