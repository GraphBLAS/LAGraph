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
    {      1, "karate.mtx" }, 
    {      1, "A.mtx" }, 
    {      1, "jagmesh7.mtx" }, 
    {      1, "ldbc-cdlp-undirected-example.mtx" }, 
    {      1, "ldbc-undirected-example.mtx" }, 
    {      1, "ldbc-wcc-example.mtx" }, 
    {      3, "LFAT5.mtx" }, 
    {      6, "LFAT5_two.mtx" }, 
    {      1, "bcsstk13.mtx" }, 
    {      1, "tree-example.mtx" }, 
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

        for (int trial = 0 ; trial <= 1 ; trial++)
        {
            // find the connected components
            OK (LAGraph_ConnectedComponents (&C, G, msg)) ;
            GxB_print (C, 2) ;

            // count the # of connected components
            int ncomponents = 0 ;
            GrB_Index n ;
            OK (GrB_Matrix_nrows (&n, G->A)) ;
            for (int i = 0 ; i < n ; i++)
            {
                int comp = -1 ;
                OK (GrB_Vector_extractElement (&comp, C, i)) ;
                if (comp == i) ncomponents++ ;
            }
            printf ("# components: %6u Matrix: %s\n", ncomponents, aname) ;
            TEST_CHECK (ncomponents == ncomp) ;

            // check the result
            OK (LG_check_cc (C, G, msg)) ;

            // convert to directed with symmetric pattern for next trial
            G->kind = LAGRAPH_ADJACENCY_DIRECTED ;
            G->A_structure_is_symmetric = LAGRAPH_TRUE ;
        }

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
