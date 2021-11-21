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
#include <LAGraphX.h>

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

//------------------------------------------------------------------------------
// count_connected_components: count the # of components in a component vector
//------------------------------------------------------------------------------

int count_connected_components (GrB_Vector C) ;

int count_connected_components (GrB_Vector C)
{
    GrB_Index n ;
    OK (GrB_Vector_size (&n, C)) ;
    int ncomponents = 0 ;
    for (int i = 0 ; i < n ; i++)
    {
        int comp = -1 ;
        int result = GrB_Vector_extractElement (&comp, C, i) ;
        if (comp == i) ncomponents++ ;
    }
    return (ncomponents) ;
}

//----------------------------------------------------------------------------
// test_cc_matrices: test with several matrices
//----------------------------------------------------------------------------

void test_cc_matrices (void)
{

#if !LG_SUITESPARSE
    printf ("SuiteSparse required for CC test\n") ;
    return ;
#endif

    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Vector C = NULL, C2 = NULL ;
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
        GrB_Index n ;
        OK (GrB_Matrix_nrows (&n, A)) ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        for (int trial = 0 ; trial <= 1 ; trial++)
        {
            // find the connected components
            OK (LAGraph_ConnectedComponents (&C, G, msg)) ;
            OK (LAGraph_Vector_print (C, 3, stdout, msg)) ;

            // count the # of connected components
            int ncomponents = count_connected_components (C) ;
            printf ("# components: %6u Matrix: %s\n", ncomponents, aname) ;
            TEST_CHECK (ncomponents == ncomp) ;

            // check the result
            OK (LG_check_cc (C, G, msg)) ;

            if (trial == 0)
            {
                // find the connected components with cc_boruvka
                printf ("BORUVKA:\n") ;
                OK (LAGraph_cc_boruvka (&C2, G->A, false)) ;
                OK (LAGraph_Vector_print (C2, 3, stdout, msg)) ;
                ncomponents = count_connected_components (C2) ;
                TEST_CHECK (ncomponents == ncomp) ;
                OK (LG_check_cc (C2, G, msg)) ;
                OK (GrB_free (&C2)) ;

                // find the connected components with cc_lacc
                #if 0
                printf ("LACC:\n") ;
                OK (LAGraph_cc_lacc (&C2, G->A, false)) ;
                // OK (LAGraph_Vector_print (C2, 3, stdout, msg)) ;
                GxB_print (C2, 3) ;
                ncomponents = count_connected_components (C2) ;
                TEST_CHECK (ncomponents == ncomp) ;
                OK (LG_check_cc (C2, G, msg)) ;
                OK (GrB_free (&C2)) ;
                printf ("did LACC\n") ;
                #endif
            }

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
