//------------------------------------------------------------------------------
// LAGraph/src/test/test_MMRead.c:  test LAGraph_MMRead and LAGraph_MMWrite
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
GrB_Type atype = NULL, btype = NULL ;
const char *name, *date ;
int ver [3] ;
GrB_Index nrows, ncols, nvals ;
#define LEN 512
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// test matrices
//------------------------------------------------------------------------------

typedef struct
{
    GrB_Index nrows ;
    GrB_Index ncols ;
    GrB_Index nvals ;
    const char *type ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    // nrows ncols nvals type         name
    {    7,    7,    30, "GrB_BOOL",  "A.mtx" },
    {    7,    7,    12, "GrB_INT32", "cover.mtx" },
    {    7,    7,    12, "GrB_BOOL",  "cover_structure.mtx" },
    { 1138, 1138,  7450, "GrB_BOOL",  "jagmesh7.mtx" },
    {    8,    8,    18, "GrB_BOOL",  "ldbc-cdlp-directed-example.mtx" },
    {    8,    8,    24, "GrB_BOOL",  "ldbc-cdlp-undirected-example.mtx" },
    {   10,   10,    17, "GrB_BOOL",  "ldbc-directed-example-bool.mtx" },
    {   10,   10,    17, "GrB_FP64",  "ldbc-directed-example.mtx" },
    {   10,   10,    17, "GrB_BOOL",  "ldbc-directed-example-unweighted.mtx" },
    {    9,    9,    24, "GrB_BOOL",  "ldbc-undirected-example-bool.mtx" },
    {    9,    9,    24, "GrB_FP64",  "ldbc-undirected-example.mtx" },
    {    9,    9,    24, "GrB_BOOL",  "ldbc-undirected-example-unweighted.mtx"},
    {   10,   10,    30, "GrB_INT64", "ldbc-wcc-example.mtx" },
    {   14,   14,    46, "GrB_FP64",  "LFAT5.mtx" },
    {    6,    6,     8, "GrB_INT64", "msf1.mtx" },
    {    8,    8,    12, "GrB_INT64", "msf2.mtx" },
    {    5,    5,     7, "GrB_INT64", "msf3.mtx" },
    {    8,    8,    28, "GrB_BOOL",  "sample2.mtx" },
    {    8,    8,    12, "GrB_BOOL",  "sample.mtx" },
    {   64,    1,    64, "GrB_INT64", "sources_7.mtx" },
    { 1000, 1000,  3996, "GrB_FP64",  "olm1000.mtx" },
    { 2003, 2003, 83883, "GrB_FP64",  "bcsstk13.mtx" },
    { 2500, 2500, 12349, "GrB_FP64",  "cryg2500.mtx" },
    {    6,    6,    10, "GrB_INT64", "tree-example.mtx" },
    {   67,   67,   294, "GrB_FP64",  "west0067.mtx" },
    {   27,   51,   102, "GrB_FP64",  "lp_afiro.mtx" },
    {   27,   51,   102, "GrB_BOOL",  "lp_afiro_structure.mtx" },
    {   34,   34,   156, "GrB_BOOL",  "karate.mtx" },
    {    7,    7,    12, "GrB_BOOL",  "matrix_bool.mtx" },
    {    7,    7,    12, "GrB_INT8",  "matrix_int8.mtx" },
    {    7,    7,    12, "GrB_INT16", "matrix_int16.mtx" },
    {    7,    7,    12, "GrB_INT32", "matrix_int32.mtx" },
    {    7,    7,    12, "GrB_INT64", "matrix_int64.mtx" },
    {    7,    7,    12, "GrB_UINT8", "matrix_uint8.mtx" },
    {    7,    7,    12, "GrB_UINT16","matrix_uint16.mtx" },
    {    7,    7,    12, "GrB_UINT32","matrix_uint32.mtx" },
    {    7,    7,    12, "GrB_UINT64","matrix_uint64.mtx" },
    {    7,    7,    12, "GrB_FP32",  "matrix_fp32.mtx" },
    {    7,    7,    12, "GrB_BOOL",  "matrix_fp32_structure.mtx" },
    {    7,    7,    12, "GrB_FP64",  "matrix_fp64.mtx" },
    {   67,   67,   294, "GrB_FP64",  "west0067_jumbled.mtx" },
    {    6,    6,    20, "GrB_FP32",  "skew_fp32.mtx" },
    {    6,    6,    20, "GrB_FP64",  "skew_fp64.mtx" },
    {    6,    6,    20, "GrB_INT8",  "skew_int8.mtx" },
    {    6,    6,    20, "GrB_INT16", "skew_int16.mtx" },
    {    6,    6,    20, "GrB_INT32", "skew_int32.mtx" },
    {    6,    6,    20, "GrB_INT64", "skew_int64.mtx" },
    {    7,    7,    12, "GrB_INT32", "structure.mtx" },
    {    3,    3,     9, "GrB_FP64",  "full.mtx" },
    {    4,    4,    16, "GrB_FP64",  "full_symmetric.mtx" },
    {    3,    4,     0, "GrB_INT32", "empty.mtx" },
    { 0, 0, 0, "", "" },
} ;

//------------------------------------------------------------------------------
// setup: start a test
//------------------------------------------------------------------------------

void setup (void)
{
    printf ("\nsetup: %s\n", __FILE__) ;
    printf ("data is in [%s]\n", LG_DATA_DIR) ;
    OK (LAGraph_Init (msg)) ;
    #if LG_SUITESPARSE
    OK (GxB_get (GxB_LIBRARY_NAME, &name)) ;
    OK (GxB_get (GxB_LIBRARY_DATE, &date)) ;
    OK (GxB_get (GxB_LIBRARY_VERSION, ver)) ;
    #endif
}

//------------------------------------------------------------------------------
// teardown: finalize a test
//------------------------------------------------------------------------------

void teardown (void)
{
    #if LG_SUITESPARSE
    printf ("\n%s %d.%d.%d (%s)\n", name, ver [0], ver [1], ver [2], date) ;
    #endif
    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    TEST_CHECK (A == NULL) ;
    TEST_CHECK (B == NULL) ;
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_MMRead:  read a set of matrices, check their stats, and write them out
//------------------------------------------------------------------------------

void test_MMRead (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth file
        //----------------------------------------------------------------------

        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\n============= %2d: %s\n", k, aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;

        //----------------------------------------------------------------------
        // check its stats
        //----------------------------------------------------------------------

        OK (GrB_Matrix_nrows (&nrows, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;
        OK (GrB_Matrix_nvals (&nvals, A)) ;
        TEST_CHECK (nrows == files [k].nrows) ;
        TEST_CHECK (ncols == files [k].ncols) ;
        TEST_CHECK (nvals == files [k].nvals) ;
        #if LG_SUITESPARSE
        OK (GxB_Matrix_type (&btype, A)) ;
        TEST_CHECK (atype == btype) ;
        #endif
        const char *tname = typename (atype) ;
        TEST_CHECK (tname != NULL) ;
        OK (strcmp (tname, files [k].type)) ;
        TEST_MSG ("Stats are wrong for %s\n", aname) ;

        //----------------------------------------------------------------------
        // pretty-print the matrix
        //----------------------------------------------------------------------

        for (int pr = 0 ; pr <= 2 ; pr++)
        {
            printf ("\nPretty-print %s: pr=%d:\n", aname, pr) ;
            OK (LAGraph_Matrix_print (A, pr, stdout, msg)) ;
        }

        //----------------------------------------------------------------------
        // write it to a temporary file
        //----------------------------------------------------------------------

        f = tmpfile ( ) ;
        OK (LAGraph_MMWrite_type (A, atype, f, NULL, msg)) ;
        TEST_MSG ("Failed to write %s to a temp file\n", aname) ;

        //----------------------------------------------------------------------
        // load it back in again
        //----------------------------------------------------------------------

        rewind (f) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        TEST_MSG ("Failed to load %s from a temp file\n", aname) ;
        OK (fclose (f)) ;       // close and delete the temporary file

        //----------------------------------------------------------------------
        // ensure A and B are the same
        //----------------------------------------------------------------------

        TEST_CHECK (atype == btype) ;
        bool ok ;
        OK (LAGraph_IsEqual_type (&ok, A, B, atype, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Failed test for equality, file: %s\n", aname) ;

        //----------------------------------------------------------------------
        // free workspace
        //----------------------------------------------------------------------

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_karate: read in karate graph from a file and compare it known graph
//-----------------------------------------------------------------------------

void test_karate (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    setup ( ) ;

    //--------------------------------------------------------------------------
    // load in the data/karate.mtx file as the matrix A
    //--------------------------------------------------------------------------

    FILE *f = fopen (LG_DATA_DIR "karate.mtx", "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
    TEST_CHECK (atype == GrB_BOOL) ;
    OK (fclose (f)) ;
    OK (LAGraph_Matrix_print_type (A, atype, 2, stdout, msg)) ;
    TEST_MSG ("Loading of A matrix failed: karate matrix") ;

    //--------------------------------------------------------------------------
    // load in the matrix defined by graph_zachary_karate.h as the matrix B
    //--------------------------------------------------------------------------

    OK (GrB_Matrix_new (&B, GrB_BOOL, ZACHARY_NUM_NODES, ZACHARY_NUM_NODES)) ;
    OK (GrB_Matrix_build (B, ZACHARY_I, ZACHARY_J, ZACHARY_V,
        ZACHARY_NUM_EDGES, GrB_LOR)) ;
    OK (LAGraph_Matrix_print_type (B, GrB_BOOL, 2, stdout, msg)) ;
    TEST_MSG ("Loading of B matrix failed: karate matrix") ;

    //--------------------------------------------------------------------------
    // ensure A and B are the same
    //--------------------------------------------------------------------------

    bool ok ;
    OK (LAGraph_IsEqual_type (&ok, A, B, GrB_BOOL, msg)) ;
    TEST_CHECK (ok) ;
    TEST_MSG ("Test for A and B equal failed: karate matrix") ;

    //--------------------------------------------------------------------------
    // free workspace and finish the test
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_failures: test for failure modes of LAGraph_MMRead and MMWrite
//-----------------------------------------------------------------------------

typedef struct
{
    int error ;
    const char *name ;
}
mangled_matrix_info ;

const mangled_matrix_info mangled_files [ ] =
{
//  error  filename              how the matrix is mangled
    -1002, "mangled1.mtx",       // bad header
    -1002, "mangled2.mtx",       // bad header
    -1002, "mangled3.mtx",       // bad type
    -1,    "complex.mtx",        // valid complex matrix, not supported
    -1002, "mangled4.mtx",       // bad format
    -1002, "mangled5.mtx",       // invalid combination of format options
    -1002, "mangled6.mtx",       // invalid combination of format options
    -1002, "mangled7.mtx",       // invalid GraphBLAS type
    -1002, "mangled8.mtx",       // invalid first line
    -1002, "mangled9.mtx",       // invalid matrix: symmetric and rectangular
    -1002, "mangled10.mtx",      // invalid matrix: truncated
    -1002, "mangled11.mtx",      // invalid matrix: entries mangled
    -1002, "mangled12.mtx",      // invalid matrix: entries mangled
    -GrB_INVALID_INDEX, "mangled13.mtx", // invalid matrix: indices out of range
    -1002, "mangled14.mtx",      // invalid matrix: duplicate entries
    -1002, "mangled_bool.mtx",   // invalid matrix: entry value out of range
    -1002, "mangled_int8.mtx",   // invalid matrix: entry value out of range
    -1002, "mangled_int16.mtx",  // invalid matrix: entry value out of range
    -1002, "mangled_int32.mtx",  // invalid matrix: entry value out of range
    -1002, "mangled_uint8.mtx",  // invalid matrix: entry value out of range
    -1002, "mangled_uint16.mtx", // invalid matrix: entry value out of range
    -1002, "mangled_uint32.mtx", // invalid matrix: entry value out of range
    -1002, "mangled_skew.mtx",   // invalid matrix: unsigned skew invalid
    0, "",
} ;

void test_MMRead_failures (void)
{
    setup ( ) ;
    printf ("\nTesting error handling of LAGraph_MMRead when giving it "
        "mangled matrices:\n") ;

    // input arguments are NULL
    TEST_CHECK (LAGraph_MMRead (NULL, NULL, NULL, msg) == -1001) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (LAGraph_MMRead (&A, NULL, NULL, msg) == -1001) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (LAGraph_MMRead (&A, &atype, NULL, msg) == -1001) ;
    printf ("msg: [%s]\n", msg) ;

    // matrix files are mangled in some way, or unsupported
    for (int k = 0 ; ; k++)
    {
        const char *aname = mangled_files [k].name ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        int error = mangled_files [k].error ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        printf ("file: [%s]\n", filename) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        int status = LAGraph_MMRead (&A, &atype, f, msg) ;
        TEST_CHECK (status == error || status == -error) ;
        if (status == error || status == -error)
        {
            printf ("    got the error we expected: %d [%s]\n", status, msg) ;
        }
        OK (fclose (f)) ;
        TEST_CHECK (A == NULL) ;
    }

    // typename is invalid
    const char *tname = typename (NULL) ;
    TEST_CHECK (tname == NULL) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_jumbled: test reading a jumbled matrix
//-----------------------------------------------------------------------------

void test_jumbled (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    setup ( ) ;

    //--------------------------------------------------------------------------
    // load in the data/west0067.mtx file as the matrix A
    //--------------------------------------------------------------------------

    FILE *f = fopen (LG_DATA_DIR "west0067.mtx", "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
    TEST_CHECK (atype == GrB_FP64) ;
    OK (fclose (f)) ;
    TEST_MSG ("Loading of west0067.mtx failed") ;

    //--------------------------------------------------------------------------
    // load in the data/west0067_jumbled.mtx file as the matrix B
    //--------------------------------------------------------------------------

    f = fopen (LG_DATA_DIR "west0067_jumbled.mtx", "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
    TEST_CHECK (btype == GrB_FP64) ;
    OK (fclose (f)) ;
    TEST_MSG ("Loading of west0067_jumbled.mtx failed") ;

    //--------------------------------------------------------------------------
    // ensure A and B are the same
    //--------------------------------------------------------------------------

    bool ok ;
    OK (LAGraph_IsEqual_type (&ok, A, B, atype, msg)) ;
    TEST_CHECK (ok) ;
    TEST_MSG ("Test for A and B equal failed: west0067_jumbled.mtx matrix") ;

    //--------------------------------------------------------------------------
    // free workspace and finish the test
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_MMWrite: test LAGraph_MMWrite
//-----------------------------------------------------------------------------

const char* files_for_MMWrite [ ] =
{
    "west0067.mtx",
    "full.mtx",
    "cover.mtx",
    ""
} ;

void test_MMWrite (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth file
        //----------------------------------------------------------------------

        const char *aname = files_for_MMWrite [k] ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\n============= %2d: %s\n", k, aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;

        //----------------------------------------------------------------------
        // create a file for comments
        //----------------------------------------------------------------------

        FILE *fcomments = fopen (LG_DATA_DIR "comments.txt", "w") ;
        TEST_CHECK (fcomments != NULL) ;
        fprintf (fcomments, " comments for %s\n", aname) ;
        fprintf (fcomments, " this file was created by test_MMRead.c\n") ;
        fclose (fcomments) ;
        TEST_MSG ("Failed to create comments.txt") ;

        //----------------------------------------------------------------------
        // write the matrix to the data/comments_*.mtx file
        //----------------------------------------------------------------------

        snprintf (filename, LEN, LG_DATA_DIR "comments_%s", aname) ;
        fcomments = fopen (LG_DATA_DIR "comments.txt", "r") ;
        FILE *foutput = fopen (filename, "w") ;
        TEST_CHECK (foutput != NULL) ;
        TEST_CHECK (fcomments != NULL) ;
        if (atype == GrB_FP64)
        {
            // select the type automatically
            OK (LAGraph_MMWrite (A, foutput, fcomments, msg)) ;
        }
        else
        {
            // pass in the type
            OK (LAGraph_MMWrite_type (A, atype, foutput, fcomments, msg)) ;
        }
        fclose (fcomments) ;
        fclose (foutput) ;
        TEST_MSG ("Failed to create %s", filename) ;

        //----------------------------------------------------------------------
        // load in the data/comments_.mtx file as the matrix B
        //----------------------------------------------------------------------

        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;
        TEST_CHECK (btype == atype) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of %s failed", filename) ;

        //----------------------------------------------------------------------
        // ensure A and B are the same
        //----------------------------------------------------------------------

        bool ok ;
        OK (LAGraph_IsEqual_type (&ok, A, B, atype, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Test for A and B equal failed: %s", filename) ;

        //----------------------------------------------------------------------
        // write a nan
        //----------------------------------------------------------------------

        if (k == 0)
        {
            OK (GrB_Matrix_setElement (A, NAN, 0, 0)) ;
            double a ;
            OK (GrB_Matrix_extractElement (&a, A, 0, 0)) ;
            TEST_CHECK (isnan (a)) ;
            foutput = fopen (filename, "w") ;
            fcomments = fopen (LG_DATA_DIR "comments.txt", "r") ;
            TEST_CHECK (foutput != NULL) ;
            OK (LAGraph_MMWrite_type (A, GrB_FP64, foutput, fcomments, msg)) ;
            fclose (fcomments) ;
            fclose (foutput) ;
            OK (GrB_free (&A)) ;
            f = fopen (filename, "r") ;
            TEST_CHECK (f != NULL) ;
            OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
            fclose (f) ;
            a = 0 ;
            OK (GrB_Matrix_extractElement (&a, A, 0, 0)) ;
            TEST_CHECK (isnan (a)) ;
        }

        //----------------------------------------------------------------------
        // free workspace
        //----------------------------------------------------------------------

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_MMWrite_failures: test error handling of LAGraph_MMWrite
//-----------------------------------------------------------------------------

typedef int mytype ;

void test_MMWrite_failures (void)
{
    setup ( ) ;
    printf ("\nTesting error handling of LAGraph_MMWrite\n") ;

    // input arguments are NULL
    TEST_CHECK (LAGraph_MMWrite (NULL, NULL, NULL, msg) == -1001) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (LAGraph_MMWrite_type (NULL, NULL, NULL, NULL, msg) == -1001) ;

    // attempt to print a matrix with a user-defined type, which should fail
    FILE *f = tmpfile ( ) ;
    TEST_CHECK (f != NULL) ;
    OK (GrB_Type_new (&atype, sizeof (mytype))) ;
    OK (GrB_Matrix_new (&A, atype, 4, 4)) ;
    int status = LAGraph_Matrix_print_type (A, atype, 3, stdout, msg) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (status == -1002) ;
    status = LAGraph_MMWrite_type (A, atype, f, NULL, msg) ;
    printf ("msg: [%s]\n", msg) ;
    TEST_CHECK (status == -1006) ;
    OK (GrB_free (&atype)) ;
    OK (GrB_free (&A)) ;
    OK (fclose (f)) ;       // close and delete the temporary file

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_MMReadWrite_brutal
//------------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_MMReadWrite_brutal (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    OK (LG_brutal_setup (msg)) ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth file
        //----------------------------------------------------------------------

        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\n============= %2d: %s\n", k, aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;
        printf ("\n") ;

        //----------------------------------------------------------------------
        // write it to a temporary file
        //----------------------------------------------------------------------

        for (int nbrutal = 0 ; ; nbrutal++)
        {
            /* allow for only nbrutal mallocs before 'failing' */
            printf (".") ;
            LG_brutal = nbrutal ;
            /* try the method with brutal malloc */
            f = tmpfile ( ) ;   // create a new temp file for each trial
            int brutal_result = LAGraph_MMWrite_type (A, atype, f, NULL, msg) ;
            if (brutal_result >= 0)
            {
                /* the method finally succeeded */
                // leave the file open for the next phase
                printf (" MMWrite ok: %d mallocs\n", nbrutal) ;
                break ;
            }
            OK (fclose (f)) ;   // close and delete the file and try again
            if (nbrutal > 10000) { printf ("Infinite!\n") ; abort ( ) ; }
        }
        LG_brutal = -1 ;  /* turn off brutal mallocs */

        //----------------------------------------------------------------------
        // load it back in again
        //----------------------------------------------------------------------

        for (int nbrutal = 0 ; ; nbrutal++)
        {
            /* allow for only nbrutal mallocs before 'failing' */
            printf (".") ;
            LG_brutal = nbrutal ;
            /* try the method with brutal malloc */
            rewind (f) ;        // rewind the temp file for each trial
            int brutal_result = LAGraph_MMRead (&B, &btype, f, msg) ;
            if (brutal_result >= 0)
            {
                /* the method finally succeeded */
                printf (" MMRead ok: %d mallocs\n", nbrutal) ;
                OK (fclose (f)) ;   // finally close and delete the temp file
                break ;
            }
            if (nbrutal > 10000) { printf ("Infinite!\n") ; abort ( ) ; }
        }
        LG_brutal = -1 ;  /* turn off brutal mallocs */

        //----------------------------------------------------------------------
        // ensure A and B are the same
        //----------------------------------------------------------------------

        TEST_CHECK (atype == btype) ;
        bool ok ;
        LG_BRUTAL (LAGraph_IsEqual_type (&ok, A, B, atype, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Failed test for equality, file: %s\n", aname) ;

        //----------------------------------------------------------------------
        // free workspace
        //----------------------------------------------------------------------

        OK (GrB_free (&A)) ;
        OK (GrB_free (&B)) ;
    }

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "MMRead", test_MMRead },
    { "karate", test_karate },
    { "MMRead_failures", test_MMRead_failures },
    { "jumbled", test_jumbled },
    { "MMWrite", test_MMWrite },
    { "MMWrite_failures", test_MMWrite_failures },
    #if LG_SUITESPARSE
    { "MMReadWrite_brutal", test_MMReadWrite_brutal },
    #endif
    { NULL, NULL }
} ;

