//----------------------------------------------------------------------------
// LAGraph/experimental/test/test_lcc.c: tests for Local Clustering Coefficient
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

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include "LG_Xtest.h"

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
void test_lcc (void)
{
    LAGraph_Init (msg) ;
    #if LAGRAPH_SUITESPARSE

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

        // construct a directed graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        // check for self-edges
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg));
        OK (LAGraph_Cached_NSelfEdges (G, msg)) ;

        GrB_Vector c = NULL ;

        // compute the local clustering coefficient
        OK (LAGraph_lcc (&c, G, msg)) ;

        GrB_Index n ;
        OK (GrB_Vector_size (&n, c)) ;
        LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT ;

        GrB_Vector cgood = NULL ;
        OK (LG_check_lcc(&cgood, G, msg)) ;
        OK (GrB_wait (cgood, GrB_MATERIALIZE)) ;
        // cgood = abs (cgood - c)
        OK (GrB_eWiseAdd (cgood, NULL, NULL, GrB_MINUS_FP64, cgood, c,
            NULL)) ;
        OK (GrB_apply (cgood, NULL, NULL, GrB_ABS_FP64, cgood, NULL)) ;
        double err = 0 ;
        // err = max (cgood)
        OK (GrB_reduce (&err, NULL, GrB_MAX_MONOID_FP64, cgood, NULL)) ;
        printf ("err: %g\n", err) ;
        TEST_CHECK (err < 1e-6) ;
        OK (GrB_free (&cgood)) ;

        printf ("\nlcc:\n") ;
        OK (LAGraph_Vector_Print (c, pr, stdout, msg)) ;
        OK (GrB_free (&c)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    #else
    printf ("test skipped\n") ;
    #endif
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors (void)
{
    LAGraph_Init (msg) ;
    #if LAGRAPH_SUITESPARSE

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    TEST_MSG ("Loading of adjacency matrix failed") ;

    // construct an undirected graph G with adjacency matrix A
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    OK (LAGraph_Cached_IsSymmetricStructure (G, msg));
    OK (LAGraph_Cached_NSelfEdges (G, msg)) ;

    GrB_Vector c = NULL ;

    // c is NULL
    int result = LAGraph_lcc (NULL, G, msg) ;
    printf ("\nresult: %d\n", result) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    OK (LAGraph_Delete (&G, msg)) ;
    #else
    printf ("test skipped\n") ;
    #endif
    LAGraph_Finalize (msg) ;
}

//****************************************************************************

TEST_LIST = {
    {"lcc", test_lcc},
    {"lcc_errors", test_errors},
    {NULL, NULL}
};
