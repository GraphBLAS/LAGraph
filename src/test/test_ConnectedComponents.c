//----------------------------------------------------------------------------
// LAGraph/src/test/test_ConnectedComponents.c: test cases for CC
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
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    uint32_t ncomponents ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    {      0, "karate.mtx" }, 
    {      0, "A.mtx" }, 
    {      0, "jagmesh7.mtx" }, 
    {      0, "ldbc-cdlp-undirected-example.mtx" }, 
    {      0, "ldbc-undirected-example.mtx" }, 
    {      0, "ldbc-wcc-example.mtx" }, 
    {      0, "LFAT5.mtx" }, 
    {      0, "LFAT5_two.mtx" }, 
    {      0, "bcsstk13.mtx" }, 
    {      0, "tree-example.mtx" }, 
    {      0, "" }, 
} ;

//----------------------------------------------------------------------------
// test_cc_matrices: test with several matrices
//----------------------------------------------------------------------------

void test_cc_matrices (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Vector C = NULL ;
    GrB_Type atype = NULL ;

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        uint32_t ncomp = files [k].ncomponents ;
        if (strlen (aname) == 0) break;
        printf ("\nMatrix: %s\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // find the connected components
        OK (LAGraph_ConnectedComponents (&C, G, msg)) ;
        GxB_print (C, 2) ;

        uint32_t ncomponents ;
        OK (GrB_reduce (&ncomponents, NULL, GrB_MAX_MONOID_UINT32, C, NULL)) ;
        printf ("# components: %6u Matrix: %s\n", ncomponents+1, aname) ;

        // TEST_CHECK (ncomponents == ncomp) ;
        int result = (LG_check_cc (C, G, msg)) ;
        printf ("msg: %s\n", msg) ;
        OK (result) ;

        OK (LAGraph_Delete (&G, msg)) ;
        OK (GrB_free (&C)) ;
    }

    LAGraph_Finalize(msg);
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"cc", test_cc_matrices},
    {NULL, NULL}
};
