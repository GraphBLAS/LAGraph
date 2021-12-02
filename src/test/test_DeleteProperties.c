//------------------------------------------------------------------------------
// LAGraph/src/test/test_DeleteProperties.c:  test LAGraph_DeleteProperties
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

LAGraph_Graph G = NULL ;
char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL ;
GrB_Type atype = NULL ;
#define LEN 512
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// setup: start a test
//------------------------------------------------------------------------------

void setup (void)
{
    OK (LAGraph_Init (msg)) ;
}

//------------------------------------------------------------------------------
// teardown: finalize a test
//------------------------------------------------------------------------------

void teardown (void)
{
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_DeleteProperties:  test LAGraph_DeleteProperties
//------------------------------------------------------------------------------

typedef struct
{
    LAGraph_Kind kind ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    LAGRAPH_ADJACENCY_DIRECTED,   "cover.mtx",
    LAGRAPH_ADJACENCY_DIRECTED,   "ldbc-directed-example.mtx",
    LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-undirected-example.mtx",
    LAGRAPH_ADJACENCY_UNDIRECTED, "A.mtx",
    LAGRAPH_ADJACENCY_UNDIRECTED, "bcsstk13.mtx",
    LAGRAPH_UNKNOWN,              ""
} ;

void test_DeleteProperties (void)
{
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        LAGraph_Kind kind = files [k].kind ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // construct the graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        TEST_CHECK (A == NULL) ;

        // create all properties (see test_Property_* for tests of content)
        OK (LAGraph_Property_RowDegree (G, msg)) ;
        OK (LAGraph_Property_ColDegree (G, msg)) ;
        OK (LAGraph_Property_AT (G, msg)) ;
        OK (LAGraph_Property_ASymmetricStructure (G, msg)) ;

        // print them
        printf ("\nGraph: ndiag %g, symmetric structure: %d\n",
            (double) G->ndiag, G->A_structure_is_symmetric) ;
        printf ("  adj matrix: ") ;
        OK (LAGraph_Matrix_print_type (G->A, atype, 2, stdout, msg)) ;
        printf ("  row degree: ") ;
        OK (LAGraph_Vector_print_type (G->rowdegree, G->rowdegree_type, 2,
            stdout, msg)) ;
        if (kind == LAGRAPH_ADJACENCY_DIRECTED)
        {
            printf ("  adj transposed: ") ;
            OK (LAGraph_Matrix_print_type (G->AT, atype, 2, stdout, msg)) ;
            printf ("  col degree: ") ;
            OK (LAGraph_Vector_print_type (G->coldegree, G->coldegree_type, 2,
                stdout, msg)) ;
        }
        else
        {
            TEST_CHECK (G->AT == NULL) ;
            TEST_CHECK (G->coldegree == NULL) ;
        }

        for (int trial = 0 ; trial <= 1 ; trial++)
        {
            // delete all the properties
            OK (LAGraph_DeleteProperties (G, msg)) ;
            TEST_CHECK (G->AT == NULL) ;
            TEST_CHECK (G->rowdegree == NULL) ;
            TEST_CHECK (G->rowdegree_type == NULL) ;
            TEST_CHECK (G->coldegree == NULL) ;
            TEST_CHECK (G->coldegree_type == NULL) ;
        }

        OK (LAGraph_Delete (&G, msg)) ;
    }

    OK (LAGraph_DeleteProperties (NULL, msg)) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_del_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        LAGraph_Kind kind = files [k].kind ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // construct the graph G with adjacency matrix A
        LG_BRUTAL (LAGraph_New (&G, &A, atype, kind, msg)) ;
        TEST_CHECK (A == NULL) ;

        // create all properties (see test_Property_* for tests of content)
        LG_BRUTAL (LAGraph_Property_RowDegree (G, msg)) ;
        LG_BRUTAL (LAGraph_Property_ColDegree (G, msg)) ;
        LG_BRUTAL (LAGraph_Property_AT (G, msg)) ;
        LG_BRUTAL (LAGraph_Property_ASymmetricStructure (G, msg)) ;

        for (int trial = 0 ; trial <= 1 ; trial++)
        {
            // delete all the properties
            LG_BRUTAL (LAGraph_DeleteProperties (G, msg)) ;
            TEST_CHECK (G->AT == NULL) ;
            TEST_CHECK (G->rowdegree == NULL) ;
            TEST_CHECK (G->rowdegree_type == NULL) ;
            TEST_CHECK (G->coldegree == NULL) ;
            TEST_CHECK (G->coldegree_type == NULL) ;
        }
    
        LG_BRUTAL (LAGraph_Delete (&G, msg)) ;
        LG_BRUTAL (LAGraph_DeleteProperties (NULL, msg)) ;
    }

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "Property_DeleteProperties", test_DeleteProperties },
    { "Property_DeleteProperties_brutal", test_del_brutal },
    { NULL, NULL }
} ;

