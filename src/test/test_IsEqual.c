//------------------------------------------------------------------------------
// LAGraph/src/test/test_IsEqual.c:  test LAGraph_IsEqual
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

int status ;
GrB_Info info ;
char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL, B = NULL ;
GrB_Type atype = NULL, btype = NULL, mytype = NULL ;
#define LEN 512
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// test matrices
//------------------------------------------------------------------------------

typedef struct
{
    bool isequal ;
    bool isequal_auto ;
    const char *typename ;
    const char *matrix1 ;
    const char *matrix2 ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    //   iseq  type          matrix1             matrix2
    {    0, 0, "GrB_BOOL"  , "A.mtx"           , "cover.mtx" },
    {    0, 0, "GrB_BOOL"  , "A.mtx"           , "A2.mtx" },
    {    1, 0, "GrB_BOOL"  , "cover.mtx"       , "cover_structure.mtx" },
    {    0, 0, "GrB_INT32" , "cover.mtx"       , "cover_structure.mtx" },
    {    1, 1, "GrB_FP64"  , "LFAT5.mtx"       , "LFAT5.mtx" },
    {    0, 0, "GrB_BOOL"  , "sample2.mtx"     , "sample.mtx" },
    {    1, 1, "GrB_BOOL"  , "sample.mtx"      , "sample.mtx" },
    {    1, 1, "GrB_FP64"  , "matrix_int32.mtx", "matrix_int32.mtx" },
    {    1, 1, "GrB_INT32" , "matrix_int32.mtx", "matrix_int32.mtx" },
    {    0, 0, "GrB_INT32" , "matrix_int32.mtx", "matrix_int64.mtx" },
    {    1, 0, "GrB_BOOL"  , "matrix_int32.mtx", "matrix_int64.mtx" },
    {    1, 1, "GrB_FP64"  , "west0067.mtx"    , "west0067_jumbled.mtx" },
    {    0, 0, "GrB_FP64"  , "LFAT5.mtx"       , "west0067.mtx" },
    {    0, 0, "GrB_FP64"  , "empty.mtx"       , "full.mtx" },
    {    0, 0, NULL        , ""                , "" }
} ;

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
// test_IsEqual: test LAGraph_IsEqual and LAGraph_IsEqual_type
//------------------------------------------------------------------------------

void test_IsEqual (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    setup ( ) ;
    printf ("\nTesting IsEqual:\n") ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth pair of files
        //----------------------------------------------------------------------

        const char *aname = files [k].matrix1 ;
        const char *bname = files [k].matrix2 ;
        const char *typename = files [k].typename ;
        const bool isequal = files [k].isequal ;
        const bool isequal_auto = files [k].isequal_auto ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("test %2d: %s %s (%s)\n", k, aname, bname, typename) ;

        // get the type
        GrB_Type type ;
        if      (strcmp (typename, "GrB_BOOL"  ) == 0) type = GrB_BOOL   ;
        else if (strcmp (typename, "GrB_INT8"  ) == 0) type = GrB_INT8   ;
        else if (strcmp (typename, "GrB_INT16" ) == 0) type = GrB_INT16  ;
        else if (strcmp (typename, "GrB_INT32" ) == 0) type = GrB_INT32  ;
        else if (strcmp (typename, "GrB_INT64" ) == 0) type = GrB_INT64  ;
        else if (strcmp (typename, "GrB_UINT8" ) == 0) type = GrB_UINT8  ;
        else if (strcmp (typename, "GrB_UINT16") == 0) type = GrB_UINT16 ;
        else if (strcmp (typename, "GrB_UINT32") == 0) type = GrB_UINT32 ;
        else if (strcmp (typename, "GrB_UINT64") == 0) type = GrB_UINT64 ;
        else if (strcmp (typename, "GrB_FP32"  ) == 0) type = GrB_FP32   ;
        else if (strcmp (typename, "GrB_FP64"  ) == 0) type = GrB_FP64   ;
        else TEST_CHECK (false) ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", bname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", bname) ;

        //----------------------------------------------------------------------
        // compare the two matrices
        //----------------------------------------------------------------------

        bool result = false ;
        OK (LAGraph_IsEqual_type (&result, A, B, type, msg)) ;
        TEST_CHECK (result == isequal) ;

        OK (LAGraph_IsEqual (&result, A, B, msg)) ;
        TEST_CHECK (result == isequal_auto) ;

        OK (LAGraph_IsEqual (&result, A, A, msg)) ;
        TEST_CHECK (result == true) ;

        OK (LAGraph_IsEqual_type (&result, A, A, type, msg)) ;
        TEST_CHECK (result == true) ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_IsEqual_failures: test error handling of LAGraph_IsEqual*
//------------------------------------------------------------------------------

typedef int myint ;

void test_IsEqual_failures (void)
{
    setup ( ) ;
    printf ("\nTest IsEqual: error handling and special cases\n") ;

    bool result = false ;
    // not a failure, but a special case:
    OK (LAGraph_IsEqual_type (&result, NULL, NULL, GrB_BOOL, msg)) ;
    TEST_CHECK (result == true) ;

    TEST_CHECK (LAGraph_IsEqual_type (NULL, NULL, NULL, NULL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_IsEqual (NULL, NULL, NULL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    OK (GrB_Matrix_new (&A, GrB_BOOL, 2, 2)) ;
    OK (GrB_Matrix_new (&B, GrB_BOOL, 2, 2)) ;

    TEST_CHECK (LAGraph_IsEqual_type (NULL, A, B, NULL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_IsEqual_type (NULL, A, B, GrB_BOOL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_IsEqual (NULL, A, B, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    OK (LAGraph_IsEqual (&result, A, B, msg)) ;
    TEST_CHECK (result == true) ;

    OK (GrB_Type_new (&mytype, sizeof (int))) ; 
    TEST_CHECK (LAGraph_IsEqual_type (&result, A, B, mytype, msg) == -1002) ;
    printf ("msg: %s\n", msg) ;

    OK (GrB_free (&mytype)) ;
    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "IsEqual", test_IsEqual },
    { "IsEqual_failures", test_IsEqual_failures },
    { NULL, NULL }
} ;

