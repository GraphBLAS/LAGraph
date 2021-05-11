//------------------------------------------------------------------------------
// LAGraph/src/test/test_Pattern.c:  test LAGraph_Pattern
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
GrB_Matrix A = NULL, B = NULL, C = NULL ;
GrB_Type atype = NULL, btype = NULL ;
#define LEN 512
char filename [LEN+1] ;

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
// test_Pattern:  test LAGraph_Pattern
//------------------------------------------------------------------------------

const char *files [ ] =
{
    "cover",
    "lp_afiro",
    "matrix_fp32",
    ""
} ;

void test_Pattern (void)
{
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        // load the valued as A
        const char *aname = files [k] ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s.mtx", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of valued matrix failed") ;

        // load the pattern as B
        snprintf (filename, LEN, LG_DATA_DIR "%s_pattern.mtx", aname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        TEST_CHECK (btype == GrB_BOOL) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of pattern matrix failed") ;

        // C = pattern (A)
        OK (LAGraph_Pattern (&C, A, msg)) ;

        // ensure B and C are the same
        bool C_and_B_are_identical ;
        OK (LAGraph_IsEqual (&C_and_B_are_identical, C, B, NULL, msg)) ;
        TEST_CHECK (C_and_B_are_identical) ;
        TEST_MSG ("Test for C and B equal failed") ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
        OK (GrB_free (&C)) ;

    }
    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Pattern_failures:  test error handling of LAGraph_Pattern
//------------------------------------------------------------------------------

void test_Pattern_failures (void)
{
    setup ( ) ;

    C = NULL ;
    TEST_CHECK (LAGraph_Pattern (NULL, NULL, msg) == -1) ;
    printf ("\nmsg: [%s]\n", msg) ;
    TEST_CHECK (LAGraph_Pattern (&C, NULL, msg) == -1) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (C == NULL) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Pattern", test_Pattern },
    { "Pattern_failures", test_Pattern_failures },
    { NULL, NULL }
} ;

