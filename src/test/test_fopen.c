//-----------------------------------------------------------------------------
// LAGraph/src/test/test_fopen.c:  test fopen
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
void test_fopen(void)
{
    char buffer [1000] ;
    TEST_MSG ("Testing fopen") ;
    FILE *f = fopen ("../data/A.mtx", "r") ;
    printf ("f %p\n", f) ;
    char *r = fgets (buffer, 512, f) ;
    TEST_CHECK (r != NULL) ;
    printf ("[%s]\n", buffer) ;
    TEST_CHECK (f != NULL) ;
    fclose (f) ;
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"fopen", test_fopen},
    {NULL, NULL}
};
