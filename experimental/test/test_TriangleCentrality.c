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

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Matrix C = NULL ;
GrB_Type atype = NULL ;
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    uint64_t ntriangles ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    {     11, "A.mtx" },
    {   2016, "jagmesh7.mtx" },
    { 342300, "bcsstk13.mtx" },
    {     45, "karate.mtx" },
    {      6, "ldbc-cdlp-undirected-example.mtx" },
    {      4, "ldbc-undirected-example-bool.mtx" },
    {      4, "ldbc-undirected-example-unweighted.mtx" },
    {      4, "ldbc-undirected-example.mtx" },
    {      5, "ldbc-wcc-example.mtx" },
    { 0, "" },
} ;

//****************************************************************************
void test_TriangleCentrality (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        uint64_t ntriangles = files [k].ntriangles ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        // C = spones (A), in FP64
        GrB_Index n ;
        OK (GrB_Matrix_nrows (&n, A)) ;
        OK (GrB_Matrix_new (&C, GrB_FP64, n, n)) ;
        OK (GrB_assign (C, A, NULL, (double) 1, GrB_ALL, n, GrB_ALL, n,
            GrB_DESC_S)) ;
        OK (GrB_free (&A)) ;
        TEST_CHECK (A == NULL) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // construct an undirected graph G with adjacency matrix A
        OK (LAGraph_New (&G, &C, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (C == NULL) ;
        // OK (GxB_Matrix_fprint (G->A, "G->A", 2, stdout)) ;

        // check for self-edges
        OK (LAGraph_Property_NDiag (G, msg)) ;
        if (G->ndiag != 0)
        {
            // remove self-edges (TODO: make this an LAGraph utility)
            printf ("graph has %ld self edges\n", G->ndiag) ;
            OK (LAGraph_DeleteDiag (G, msg)) ;
            printf ("now has %ld self edges\n", G->ndiag) ;
            TEST_CHECK (G->ndiag == 0) ;
        }

        // compute the triangle centrality
        GrB_Vector c = NULL ;
        uint64_t ntri ;
        int retval = LAGraph_VertexCentrality_Triangle (&c, &ntri, G, msg) ;
        TEST_CHECK (retval == 0) ;
        // TEST_MSG ("retval = %d (%s)", retval, msg) ;
        printf ("# of triangles: %lu\n", ntri) ;
        TEST_CHECK (ntri == ntriangles) ;

        #if LG_SUITESPARSE
        int pr = (n <= 100) ? GxB_COMPLETE : GxB_SHORT ;
        OK (GxB_Vector_fprint (c, "centrality", pr, stdout)) ;
        #endif
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
