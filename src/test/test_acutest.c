//-----------------------------------------------------------------------------
// LAGraph/src/test/test_mmread.cpp:  test cases for LAGraph_mmread()
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include <LAGraph.h>

#include <acutest.h>

//****************************************************************************
void test_dummy(void)
{
    TEST_MSG("Testing equality %d", 42);
    TEST_CHECK(42 == 42);
    //BOOST_CHECK_EQUAL(42, 42);
}

//****************************************************************************
void test_dummy_fails(void)
{
    TEST_CHECK(42 == 0);
    //BOOST_CHECK_EQUAL(42, 0);
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"dummy", test_dummy},
    {"dummy_fails", test_dummy_fails},
    {NULL, NULL}
};
