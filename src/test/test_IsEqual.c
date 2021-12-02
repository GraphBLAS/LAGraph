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
GrB_Vector u = NULL, v = NULL ;
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
    bool isequal0 ;
    bool isequal0_auto ;
    const char *typename ;
    const char *matrix1 ;
    const char *matrix2 ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    //   iseq        type          matrix1             matrix2
    {    0, 0, 0, 0, "GrB_BOOL"  , "A.mtx"           , "cover.mtx" },
    {    0, 0, 1, 1, "GrB_BOOL"  , "A.mtx"           , "A2.mtx" },
    {    1, 0, 1, 0, "GrB_BOOL"  , "cover.mtx"       , "cover_structure.mtx" },
    {    0, 0, 0, 0, "GrB_INT32" , "cover.mtx"       , "cover_structure.mtx" },
    {    1, 1, 1, 1, "GrB_FP64"  , "LFAT5.mtx"       , "LFAT5.mtx" },
    {    0, 0, 0, 0, "GrB_BOOL"  , "sample2.mtx"     , "sample.mtx" },
    {    1, 1, 1, 1, "GrB_BOOL"  , "sample.mtx"      , "sample.mtx" },
    {    1, 1, 1, 1, "GrB_FP64"  , "matrix_int32.mtx", "matrix_int32.mtx" },
    {    1, 1, 1, 1, "GrB_INT32" , "matrix_int32.mtx", "matrix_int32.mtx" },
    {    0, 0, 0, 0, "GrB_INT32" , "matrix_int32.mtx", "matrix_int64.mtx" },
    {    1, 0, 1, 0, "GrB_BOOL"  , "matrix_int32.mtx", "matrix_int64.mtx" },
    {    1, 1, 1, 1, "GrB_FP64"  , "west0067.mtx"    , "west0067_jumbled.mtx" },
    {    1, 1, 1, 1, "GrB_FP64"  , "west0067.mtx"    , "west0067_noheader.mtx"},
    {    0, 0, 0, 0, "GrB_FP64"  , "LFAT5.mtx"       , "west0067.mtx" },
    {    0, 0, 0, 0, "GrB_FP64"  , "empty.mtx"       , "full.mtx" },
    {    1, 1, 1, 1, "GrB_FP64"  , "full.mtx"        , "full_noheader.mtx" },
    {    0, 0, 0, 0, NULL        , ""                , "" }
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
        const bool isequal0 = files [k].isequal0 ;
        const bool isequal0_auto = files [k].isequal0_auto ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("test %2d: %s %s (%s)\n", k, aname, bname, typename) ;

        // get the type
        GrB_Type type = NULL ;
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

        TEST_CHECK (type != NULL) ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;
        GrB_Index ancols ;
        OK (GrB_Matrix_ncols (&ancols, A)) ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", bname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", bname) ;
        GrB_Index bncols ;
        OK (GrB_Matrix_ncols (&bncols, B)) ;

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

        //----------------------------------------------------------------------
        // compare two vectors
        //----------------------------------------------------------------------

        OK (GrB_Vector_new (&u, atype, ancols)) ;
        OK (GrB_Vector_new (&v, atype, bncols)) ;
        OK (GrB_Col_extract (u, NULL, NULL, A, GrB_ALL, ancols, 0,
            GrB_DESC_T0)) ;
        OK (GrB_Col_extract (v, NULL, NULL, B, GrB_ALL, bncols, 0,
            GrB_DESC_T0)) ;

        OK (LAGraph_Vector_IsEqual_type (&result, u, v, type, msg)) ;
        TEST_CHECK (result == isequal0) ;

        OK (LAGraph_Vector_IsEqual (&result, u, v, msg)) ;
        TEST_CHECK (result == isequal0_auto) ;

        OK (LAGraph_Vector_IsEqual (&result, u, u, msg)) ;
        TEST_CHECK (result == true) ;

        OK (LAGraph_Vector_IsEqual_type (&result, u, u, type, msg)) ;
        TEST_CHECK (result == true) ;

        OK (GrB_free (&u)) ;
        OK (GrB_free (&v)) ;
        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_IsEqual_brutal:
//------------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_IsEqual_brutal (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    OK (LG_brutal_setup (msg)) ;
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
        const bool isequal0 = files [k].isequal0 ;
        const bool isequal0_auto = files [k].isequal0_auto ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("test %2d: %s %s (%s)\n", k, aname, bname, typename) ;

        // get the type
        GrB_Type type = NULL ;
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

        TEST_CHECK (type != NULL) ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;
        GrB_Index ancols ;
        OK (GrB_Matrix_ncols (&ancols, A)) ;

        snprintf (filename, LEN, LG_DATA_DIR "%s", bname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", bname) ;
        GrB_Index bncols ;
        OK (GrB_Matrix_ncols (&bncols, B)) ;

        //----------------------------------------------------------------------
        // compare the two matrices
        //----------------------------------------------------------------------

        bool result = false ;
        LG_BRUTAL (LAGraph_IsEqual_type (&result, A, B, type, msg)) ;
        TEST_CHECK (result == isequal) ;

        LG_BRUTAL (LAGraph_IsEqual (&result, A, B, msg)) ;
        TEST_CHECK (result == isequal_auto) ;

        LG_BRUTAL (LAGraph_IsEqual (&result, A, A, msg)) ;
        TEST_CHECK (result == true) ;

        LG_BRUTAL (LAGraph_IsEqual_type (&result, A, A, type, msg)) ;
        TEST_CHECK (result == true) ;

        //----------------------------------------------------------------------
        // compare two vectors
        //----------------------------------------------------------------------

        LG_BRUTAL (GrB_Vector_new (&u, atype, ancols)) ;
        LG_BRUTAL (GrB_Vector_new (&v, atype, bncols)) ;
        LG_BRUTAL (GrB_Col_extract (u, NULL, NULL, A, GrB_ALL, ancols, 0,
            GrB_DESC_T0)) ;
        LG_BRUTAL (GrB_Col_extract (v, NULL, NULL, B, GrB_ALL, bncols, 0,
            GrB_DESC_T0)) ;

        LG_BRUTAL (LAGraph_Vector_IsEqual_type (&result, u, v, type, msg)) ;
        TEST_CHECK (result == isequal0) ;

        LG_BRUTAL (LAGraph_Vector_IsEqual (&result, u, v, msg)) ;
        TEST_CHECK (result == isequal0_auto) ;

        LG_BRUTAL (LAGraph_Vector_IsEqual (&result, u, u, msg)) ;
        TEST_CHECK (result == true) ;

        LG_BRUTAL (LAGraph_Vector_IsEqual_type (&result, u, u, type, msg)) ;
        TEST_CHECK (result == true) ;

        LG_BRUTAL (GrB_free (&u)) ;
        OK (GrB_free (&v)) ;
        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    OK (LG_brutal_teardown (msg)) ;
}
#endif

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

    OK (LAGraph_Vector_IsEqual_type (&result, NULL, NULL, GrB_BOOL, msg)) ;
    TEST_CHECK (result == true) ;

    TEST_CHECK (LAGraph_IsEqual_type (NULL, NULL, NULL, NULL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_IsEqual (NULL, NULL, NULL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    OK (GrB_Matrix_new (&A, GrB_BOOL, 2, 2)) ;
    OK (GrB_Matrix_new (&B, GrB_BOOL, 2, 2)) ;

    OK (GrB_Vector_new (&u, GrB_BOOL, 2)) ;
    OK (GrB_Vector_new (&v, GrB_BOOL, 2)) ;

    TEST_CHECK (LAGraph_IsEqual_type (NULL, A, B, NULL, msg) == -1001) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_Vector_IsEqual_type (NULL, u, v, NULL, msg) == -1001) ;
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

    TEST_CHECK (LAGraph_Vector_IsEqual_type (&result, u, v, mytype, msg)
        == -1002) ;
    printf ("msg: %s\n", msg) ;

    OK (GrB_free (&mytype)) ;
    OK (GrB_free (&u)) ;
    OK (GrB_free (&v)) ;
    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Vector_IsEqual: test LAGraph_Vector_isEqual
//------------------------------------------------------------------------------

void test_Vector_IsEqual (void)
{
    setup ( ) ;

    bool result = false ;
    OK (LAGraph_Vector_IsEqual_op (&result, NULL, NULL, GrB_EQ_BOOL, msg)) ;
    TEST_CHECK (result == true) ;

    OK (GrB_Vector_new (&u, GrB_BOOL, 3)) ;
    OK (GrB_Vector_new (&v, GrB_BOOL, 2)) ;

    OK (LAGraph_Vector_IsEqual_op (&result, u, v, GrB_EQ_BOOL, msg)) ;
    TEST_CHECK (result == false) ;

    OK (GrB_free (&u)) ;
    OK (GrB_Vector_new (&u, GrB_BOOL, 2)) ;

    OK (LAGraph_Vector_IsEqual_op (&result, u, v, GrB_EQ_BOOL, msg)) ;
    TEST_CHECK (result == true) ;

    OK (GrB_Vector_setElement (u, true, 0)) ;
    OK (GrB_Vector_setElement (v, true, 1)) ;
    OK (LAGraph_Vector_IsEqual_op (&result, u, v, GrB_EQ_BOOL, msg)) ;
    TEST_CHECK (result == false) ;

    OK (LAGraph_Vector_IsEqual_type (&result, u, v, GrB_BOOL, msg)) ;
    TEST_CHECK (result == false) ;

    OK (GrB_free (&u)) ;
    OK (GrB_free (&v)) ;

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//------------------------------------------------------------------------------

TEST_LIST =
{
    { "IsEqual", test_IsEqual },
    { "Vector_IsEqual", test_Vector_IsEqual },
    { "IsEqual_failures", test_IsEqual_failures },
    { "IsEqual_brutal", test_IsEqual_brutal },
    { NULL, NULL }
} ;

