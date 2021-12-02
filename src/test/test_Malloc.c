//-----------------------------------------------------------------------------
// LAGraph/src/test/test_Malloc.c: test LAGraph_Malloc and related methods
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include "LAGraph_test.h"

//-----------------------------------------------------------------------------
// test_malloc
//-----------------------------------------------------------------------------

void test_malloc (void)
{
    char msg [LAGRAPH_MSG_LEN] ;
    OK (LAGraph_Init (msg)) ;

    char *p = LAGraph_Malloc (42, sizeof (char)) ;
    TEST_CHECK (p != NULL) ;
    for (int k = 0 ; k < 42 ; k++)
    {
        p [k] = (char) k ;
    }
    LAGraph_Free ((void **) &p) ;
    TEST_CHECK (p == NULL) ;

    p = LAGraph_Malloc (GrB_INDEX_MAX + 1, sizeof (char)) ;
    TEST_CHECK (p == NULL) ;

    p = LAGraph_Calloc (GrB_INDEX_MAX + 1, sizeof (char)) ;
    TEST_CHECK (p == NULL) ;

    p = LAGraph_Calloc (42, sizeof (char)) ;
    for (int k = 0 ; k < 42 ; k++)
    {
        TEST_CHECK (*p == '\0') ;
    }
    LAGraph_Free ((void **) &p) ;
    TEST_CHECK (p == NULL) ;

    LAGraph_Calloc_function = NULL ;

    p = LAGraph_Calloc (42, sizeof (char)) ;
    TEST_CHECK (p != NULL) ;
    for (int k = 0 ; k < 42 ; k++)
    {
        TEST_CHECK (*p == '\0') ;
    }

    bool ok = false ;
    char *pnew = LAGraph_Realloc (100, 42, sizeof (char), p, &ok) ;
    p = NULL ;
    TEST_CHECK (ok) ;
    TEST_CHECK (pnew != NULL) ;
    for (int k = 0 ; k < 42 ; k++)
    {
        TEST_CHECK (*pnew == '\0') ;
    }
    for (int k = 42 ; k < 100 ; k++)
    {
        pnew [k] = (char) k ;
    }
    LAGraph_Free ((void **) &pnew) ;
    TEST_CHECK (pnew == NULL) ;

    p = LAGraph_Realloc (80, 0, sizeof (char), NULL, &ok) ;
    TEST_CHECK (ok) ;
    TEST_CHECK (p != NULL) ;
    for (int k = 0 ; k < 80 ; k++)
    {
        p [k] = (char) k ;
    }

    pnew = LAGraph_Realloc (GrB_INDEX_MAX+1, 80, sizeof (char), p, &ok) ;
    TEST_CHECK (pnew == NULL) ;
    TEST_CHECK (!ok) ;

    pnew = LAGraph_Realloc (80, 80, sizeof (char), p, &ok) ;
    TEST_CHECK (pnew == p) ;
    TEST_CHECK (ok) ;
    for (int k = 0 ; k < 80 ; k++)
    {
        TEST_CHECK (pnew [k] == (char) k) ;
    }

    LAGraph_Realloc_function = NULL ;

    pnew = LAGraph_Realloc (100, 80, sizeof (char), p, &ok) ;
    TEST_CHECK (ok) ;
    for (int k = 0 ; k < 80 ; k++)
    {
        TEST_CHECK (pnew [k] == (char) k) ;
    }

    LAGraph_Free ((void **) &pnew) ;
    TEST_CHECK (pnew == NULL) ;

    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST = {
    {"test_malloc", test_malloc},
    // no brutal test needed
    {NULL, NULL}
};
