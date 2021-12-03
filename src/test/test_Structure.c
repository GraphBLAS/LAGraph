//------------------------------------------------------------------------------
// LAGraph/src/test/test_Structure.c:  test LAGraph_Structure
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
// test_Structure:  test LAGraph_Structure
//------------------------------------------------------------------------------

const char *files [ ] =
{
    "cover",
    "lp_afiro",
    "matrix_fp32",
    ""
} ;

void test_Structure (void)
{
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {
        // load the valued matrix as A
        const char *aname = files [k] ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s.mtx", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of valued matrix failed") ;

        // load the structure as B
        snprintf (filename, LEN, LG_DATA_DIR "%s_structure.mtx", aname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        TEST_CHECK (btype == GrB_BOOL) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of structure matrix failed") ;

        // C = structure (A)
        OK (LAGraph_Structure (&C, A, msg)) ;

        // ensure B and C are the same
        bool ok ;
        OK (LAGraph_IsEqual_type (&ok, C, B, GrB_BOOL, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Test for C and B equal failed") ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
        OK (GrB_free (&C)) ;
    }
    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Structure_brutal
//------------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_Structure_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    for (int k = 0 ; ; k++)
    {
        // load the valued matrix as A
        const char *aname = files [k] ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s.mtx", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of valued matrix failed") ;

        // load the structure as B
        snprintf (filename, LEN, LG_DATA_DIR "%s_structure.mtx", aname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        TEST_CHECK (btype == GrB_BOOL) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of structure matrix failed") ;

        // C = structure (A)
        LG_BRUTAL (LAGraph_Structure (&C, A, msg)) ;

        // ensure B and C are the same
        bool ok ;
        OK (LAGraph_IsEqual_type (&ok, C, B, GrB_BOOL, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Test for C and B equal failed") ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
        OK (GrB_free (&C)) ;
    }
    OK (LG_brutal_teardown (msg)) ;
}
#endif

//------------------------------------------------------------------------------
// test_Structure_failures:  test error handling of LAGraph_Structure
//------------------------------------------------------------------------------

void test_Structure_failures (void)
{
    setup ( ) ;

    C = NULL ;
    TEST_CHECK (LAGraph_Structure (NULL, NULL, msg) == -1) ;
    printf ("\nmsg: [%s]\n", msg) ;
    TEST_CHECK (LAGraph_Structure (&C, NULL, msg) == -1) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (C == NULL) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Structure", test_Structure },
    { "Structure_failures", test_Structure_failures },
    #if LG_SUITESPARSE
    { "Structure_brutal", test_Structure_brutal },
    #endif
    { NULL, NULL }
} ;

