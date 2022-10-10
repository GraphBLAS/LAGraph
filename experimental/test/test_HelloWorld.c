//----------------------------------------------------------------------------
// LAGraph/src/test/test_HelloWorld.c: test cases for LAGraph_HelloWorld
//----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

// This is a very simple "hello world" example of a test program for an
// algorithm in the experimental/algorithm folder.

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;

#define LEN 512
char filename [LEN+1] ;

void test_HelloWorld (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Init (msg) ;
    GrB_Matrix Y = NULL, A = NULL ;

    //--------------------------------------------------------------------------
    // test with the west0067 matrix
    //--------------------------------------------------------------------------

    // create the graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "west0067.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A

    // test the algorithm
    OK (LAGraph_HelloWorld (&Y, G, msg)) ;

    // print the result
    printf ("\nOutput of LAGraph_HelloWorld:\n") ;
    OK (LAGraph_Matrix_Print (Y, LAGraph_COMPLETE, stdout, msg)) ;

    // check the result (ensure Y is equal to G->A)
    bool ok ;
    OK (LAGraph_Matrix_IsEqual (&ok, Y, G->A, msg)) ;
    TEST_CHECK (ok) ;

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&Y)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    LAGraph_Finalize (msg) ;
}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"HelloWorld", test_HelloWorld},    // just one test in this example
    {NULL, NULL}
} ;

