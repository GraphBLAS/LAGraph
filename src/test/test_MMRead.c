//-----------------------------------------------------------------------------
// LAGraph/src/test/test_MMRead.c:  test cases for LAGraph_MMRead
//-----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include <LAGraph.h>
#include <acutest.h>

#define OK(method) TEST_CHECK (method == 0)

int status ;
GrB_Info info ;
char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Type atype = NULL ;
int version [3] ;
const char *date, *name ;
GrB_Index nrows, ncols, nvals ;

#define DATA_DIR "../data/"

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
    {    9,    9,    24, "GrB_FP64",  "ldbc-undirected-example-unweighted.mtx"},
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

#define LEN 512
char filename [LEN+1] ;

void setup (void)
{
    printf ("\nsetup: %s\n", __FILE__) ;
    printf ("data is in [%s]\n", DATA_DIR) ;
    OK (LAGraph_Init (msg)) ;
    OK (GxB_get (GxB_LIBRARY_NAME, &name)) ;
    OK (GxB_get (GxB_LIBRARY_DATE, &date)) ;
    OK (GxB_get (GxB_LIBRARY_VERSION, version)) ;
    printf ("%s %d.%d.%d (%s)\n", name,
        version [0], version [1], version [2], date) ;
}

void teardown (void)
{
    printf ("%s %d.%d.%d (%s)\n", name,
        version [0], version [1], version [2], date) ;
    OK (LAGraph_Delete (&G, msg)) ;
    TEST_CHECK (G == NULL) ;
    OK (GrB_free (&A)) ;
    TEST_CHECK (A == NULL) ;
    LAGraph_Finalize (msg) ;
    printf ("bye\n") ;
}

void test_MMRead (void)
{
    FILE *f = NULL ;
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        printf ("\n============= %2d: %s\n", k, aname) ;
        snprintf (filename, LEN, DATA_DIR "%s", aname) ;
        f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        OK (GrB_wait (&A)) ;
        OK (GrB_Matrix_nrows (&nrows, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;
        OK (GrB_Matrix_nvals (&nvals, A)) ;
        TEST_CHECK (nrows == files [k].nrows) ;
        TEST_CHECK (ncols == files [k].ncols) ;
        TEST_CHECK (nvals == files [k].nvals) ;
        OK (GxB_print (A, 2)) ;
        OK (GrB_free (&A)) ;
    }

    teardown ( ) ;
}

//****************************************************************************
//****************************************************************************

TEST_LIST =
{
    { "MMRead", test_MMRead },
    { NULL, NULL }
} ;

