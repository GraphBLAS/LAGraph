//----------------------------------------------------------------------------
// LAGraph/src/test/test_cdlp.c: test cases for CDLP
// ----------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//-----------------------------------------------------------------------------

// todo: write a simple cdlp method, as LG_check_cdlp, and compare its results
// with LAGraph_cdlp

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    { "A.mtx" },
    { "jagmesh7.mtx" },
    { "west0067.mtx" }, // unsymmetric
    { "bcsstk13.mtx" },
    { "karate.mtx" },
    { "ldbc-cdlp-undirected-example.mtx" },
    { "ldbc-undirected-example-bool.mtx" },
    { "ldbc-undirected-example-unweighted.mtx" },
    { "ldbc-undirected-example.mtx" },
    { "ldbc-wcc-example.mtx" },
    { "" },
} ;

//****************************************************************************
void test_cdlp (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; k++)
    {
        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;

        // construct a directed graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        // check for self-edges
        OK (LAGraph_DeleteSelfEdges (G, msg)) ;
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;

        GrB_Vector c = NULL ;

        // compute the communities with LAGraph_cdlp
        OK (LAGraph_cdlp (&c, G, 100, msg)) ;
        GrB_Index n ;
        OK (GrB_Vector_size (&n, c)) ;
        LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT ;
        printf ("\ncdlp (computed result):\n") ;
        OK (LAGraph_Vector_Print (c, pr, stdout, msg)) ;

        // compute with another method
        GrB_Vector cgood = NULL ;
        OK (LAGraph_cdlp_withsort(&cgood, G, 100, msg)) ;
        OK (GrB_wait (cgood, GrB_MATERIALIZE)) ;
        printf ("\ncdlp (known result):\n") ;
        OK (LAGraph_Vector_Print (cgood, pr, stdout, msg)) ;
        bool ok = false ;
        OK (LAGraph_Vector_IsEqual (&ok, c, cgood, msg)) ;
        TEST_CHECK (ok) ;
        OK (GrB_free (&cgood)) ;

        printf ("\ncdlp:\n") ;
        OK (LAGraph_Vector_Print (c, pr, stdout, msg)) ;
        OK (GrB_free (&c)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors (void)
{
    LAGraph_Init (msg) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    TEST_MSG ("Loading of adjacency matrix failed") ;
    fclose (f) ;

    // construct an undirected graph G with adjacency matrix A
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    OK (LAGraph_DeleteSelfEdges (G, msg)) ;
    OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;

    GrB_Vector c = NULL ;

    // c is NULL
    int result = LAGraph_cdlp (NULL, G, 100, msg) ;
    printf ("\nresult: %d\n", result) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    OK (LAGraph_Delete (&G, msg)) ;
    LAGraph_Finalize (msg) ;
}

//****************************************************************************

TEST_LIST = {
    {"cdlp", test_cdlp},
    {"cdlp_errors", test_errors},
    {NULL, NULL}
};
