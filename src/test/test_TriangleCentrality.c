//----------------------------------------------------------------------------
// LAGraph/src/test/test_TriangleCentrality.cpp: test cases for triangle
// centrality
// ----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>

#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Type atype = NULL ;
GxB_Scalar thunk = NULL ;
#define LEN 512
char filename [LEN+1] ;

const char *files [ ] =
{
    "A.mtx",
    "jagmesh7.mtx",
    "bcsstk13.mtx",
    "karate.mtx",
    "ldbc-cdlp-undirected-example.mtx",
    "ldbc-undirected-example-bool.mtx",
    "ldbc-undirected-example-unweighted.mtx",
    "ldbc-undirected-example.mtx",
    "ldbc-wcc-example.mtx",
    ""
} ;

//****************************************************************************
void test_TriangleCentrality (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k] ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // construct an undirected graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        // check for self-edges
        OK (LAGraph_Property_NDiag (G, msg)) ;
        if (G->ndiag != 0)
        {
            // remove self-edges (TODO: make this an LAGraph utility)
            printf ("graph has %ld self edges\n", G->ndiag) ;
            OK (GxB_Scalar_new (&thunk, GrB_INT64)) ;
            OK (GxB_Scalar_setElement (thunk, 0)) ;
            OK (GxB_select (G->A, NULL, NULL, GxB_OFFDIAG, G->A, thunk, NULL)) ;
            OK (GrB_free (&thunk)) ;
            // OK (LAGraph_DisplayGraph (G, 3, stdout, msg)) ;
            G->ndiag = LAGRAPH_UNKNOWN ;
            OK (LAGraph_Property_NDiag (G, msg)) ;
            printf ("now has %ld self edges\n", G->ndiag) ;
            TEST_CHECK (G->ndiag == 0) ;
        }

        // compute the triangle centrality
        GrB_Vector c = NULL ;

        int retval = LAGraph_VertexCentrality_Triangle (&c, G, msg) ;
        TEST_CHECK (retval == 0) ;
        TEST_MSG ("retval = %d (%s)", retval, msg) ;

        GrB_Index n ;
        OK (GrB_Matrix_nrows (&n, G->A)) ;
        int pr = (n <= 100) ? GxB_COMPLETE : GxB_SHORT ;
        OK (GxB_Vector_fprint (c, "centrality", pr, stdout)) ;
        OK (GrB_free (&c)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize (msg) ;
}

//****************************************************************************

TEST_LIST = {
    {"TriangleCentrality", test_TriangleCentrality},
    {NULL, NULL}
};
