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
int ver [3] ;
const char *date, *name ;
GrB_Index nrows, ncols, nvals ;
#define LEN 512
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// test matrices
//------------------------------------------------------------------------------

#define TEMP_DIR "/tmp/"

typedef struct
{
    GrB_Index nrows ;
    GrB_Index ncols ;
    GrB_Index nvals ;
    const char *type ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = {
//    nrows ncols  nvals type         name
    {    7,    7,    30, "GrB_BOOL",  "A.mtx" },
    {    7,    7,    12, "GrB_INT32", "cover.mtx" },
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
    {   34,   34,   156, "GrB_BOOL",  "karate.mtx" },
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
    OK (GxB_get (GxB_LIBRARY_NAME, &name)) ;
    OK (GxB_get (GxB_LIBRARY_DATE, &date)) ;
    OK (GxB_get (GxB_LIBRARY_VERSION, ver)) ;
    printf ("%s %d.%d.%d (%s)\n", name, ver [0], ver [1], ver [2], date) ;
}

//------------------------------------------------------------------------------
// teardown: finalize a test
//------------------------------------------------------------------------------

void teardown (void)
{
    printf ("%s %d.%d.%d (%s)\n", name, ver [0], ver [1], ver [2], date) ;
    OK (GrB_free (&A)) ;
    TEST_CHECK (A == NULL) ;
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
        printf ("\n============= %2d: %s\n", k, aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;

        //----------------------------------------------------------------------
        // check its stats
        //----------------------------------------------------------------------

        OK (GrB_Matrix_nrows (&nrows, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;
        OK (GrB_Matrix_nvals (&nvals, A)) ;
        TEST_CHECK (nrows == files [k].nrows) ;
        TEST_CHECK (ncols == files [k].ncols) ;
        TEST_CHECK (nvals == files [k].nvals) ;
        OK (GxB_Matrix_type (&btype, A)) ;
        TEST_CHECK (atype == btype) ;
        OK (GxB_print (A, 2)) ;

        const char *tname = typename (atype) ;
        TEST_CHECK (tname != NULL) ;
        OK (strcmp (tname, files [k].type)) ;

        //----------------------------------------------------------------------
        // write it to a temporary file
        //----------------------------------------------------------------------

        f = tmpfile ( ) ;
        OK (LAGraph_MMWrite (A, f, msg)) ;

        //----------------------------------------------------------------------
        // load it back in again
        //----------------------------------------------------------------------

        rewind (f) ;
        OK (LAGraph_MMRead (&B, &btype, f, msg)) ;

        //----------------------------------------------------------------------
        // close and delete the temporary file
        //----------------------------------------------------------------------

        OK (fclose (f)) ;

        //----------------------------------------------------------------------
        // ensure A and B are the same
        //----------------------------------------------------------------------

        TEST_CHECK (atype == btype) ;
        bool A_and_B_are_identical ;
        OK (LAGraph_IsEqual (&A_and_B_are_identical, A, B, NULL, msg)) ;
        TEST_CHECK (A_and_B_are_identical) ;

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
    OK (GxB_print (A, 2)) ;

    //--------------------------------------------------------------------------
    // load in the matrix defined by graph_zachary_karate.h as the matrix B
    //--------------------------------------------------------------------------

    OK (GrB_Matrix_new (&B, GrB_BOOL, ZACHARY_NUM_NODES, ZACHARY_NUM_NODES)) ;
    OK (GrB_Matrix_build (B, ZACHARY_I, ZACHARY_J, ZACHARY_V,
        ZACHARY_NUM_EDGES, GrB_LOR)) ;
    OK (GxB_print (B, 2)) ;

    //--------------------------------------------------------------------------
    // ensure A and B are the same
    //--------------------------------------------------------------------------

    bool A_and_B_are_identical ;
    OK (LAGraph_IsEqual (&A_and_B_are_identical, A, B, NULL, msg)) ;
    TEST_CHECK (A_and_B_are_identical) ;

    //--------------------------------------------------------------------------
    // free workspace and finish the test
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    OK (GrB_free (&B)) ;
    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "MMRead", test_MMRead },
    { "karate", test_karate },
    { NULL, NULL }
} ;

