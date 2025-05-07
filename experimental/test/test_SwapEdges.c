//----------------------------------------------------------------------------
// LAGraph/src/test/test_SwapEdges.c: test cases for LAGraph_HelloWorld
//----------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

//-----------------------------------------------------------------------------



#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LG_internal.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL, G_new = NULL;
GrB_Matrix A = NULL, C = NULL, A_new = NULL, C_new = NULL; 
#define LEN 512
char filename [LEN+1] ;

const char* tests [ ] =
{
    "random_unweighted_general1.mtx",
    "random_unweighted_general2.mtx",
    "random_weighted_general1.mtx",
    "random_weighted_general2.mtx",
    "bcsstk13.mtx",
    "test_FW_2500.mtx",
    ""
} ;
void test_SwapEdges (void)
{
    #if USING_GRAPHBLAS_V10
    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------
    OK (LAGraph_Init (msg)) ;

    for (int k = 0 ; ; k++)
    {
        //The following code taken from MIS tester
        // load the matrix as A
        const char *aname = tests [k];
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of valued matrix failed") ;
        printf ("\nMatrix: %s\n", aname) ;

        // C = structure of A
        OK (LAGraph_Matrix_Structure (&C, A, msg)) ;
        OK (GrB_free (&A)) ;

        // construct a directed graph G with adjacency matrix C
        OK (LAGraph_New (&G, &C, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (C == NULL) ;

        // check if the pattern is symmetric
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;

        if (G->is_symmetric_structure == LAGraph_FALSE)
        {
            // make the adjacency matrix symmetric
            OK (LAGraph_Cached_AT (G, msg)) ;
            OK (GrB_eWiseAdd (G->A, NULL, NULL, GrB_LOR, G->A, G->AT, NULL)) ;
            G->is_symmetric_structure = LAGraph_TRUE ;
        }
        G->kind = LAGraph_ADJACENCY_UNDIRECTED ;

        // check for self-edges
        OK (LAGraph_Cached_NSelfEdges (G, msg)) ;
        if (G->nself_edges != 0)
        {
            // remove self-edges
            printf ("graph has %g self edges\n", (double) G->nself_edges) ;
            OK (LAGraph_DeleteSelfEdges (G, msg)) ;
            printf ("now has %g self edges\n", (double) G->nself_edges) ;
            TEST_CHECK (G->nself_edges == 0) ;
        }

        // compute the row degree
        GrB_Index n = 0;
        OK (LAGraph_Cached_OutDegree (G, msg)) ;
        OK (GrB_Matrix_nrows(&n, G->A));

        //----------------------------------------------------------------------
        // test the algorithm
        //----------------------------------------------------------------------

        GrB_set (GrB_GLOBAL, (int32_t) (true), GxB_BURBLE) ;
        OK(LAGraph_SwapEdges( &A_new, G, (GrB_Index) 100, msg));
        GrB_set (GrB_GLOBAL, (int32_t) (false), GxB_BURBLE) ;
        printf ("Test ends:\n") ;
        printf ("%s\n", msg + 1) ;

        //----------------------------------------------------------------------
        // check results
        //----------------------------------------------------------------------
        bool ok = false;
        //Make sure we got a symetric back out:
        OK (LAGraph_New (&G_new, &A_new, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        OK (LAGraph_Cached_AT (G_new, msg)) ;
        OK (LAGraph_Matrix_IsEqual (&ok, G_new->AT, G_new->A, msg)) ;
        TEST_CHECK (ok) ;

        //Make sure no self edges created.
        OK (LAGraph_Cached_NSelfEdges (G_new, msg)) ;
        TEST_CHECK (G_new->nself_edges == 0);

        // Check nvals stay the same. 
        GrB_Index edge_count, new_edge_count;
        OK (GrB_Matrix_nvals(&edge_count, G->A)) ;
        OK (GrB_Matrix_nvals(&new_edge_count, G_new->A)) ;
        printf("old: %ld, new: %ld", edge_count,new_edge_count);
        TEST_CHECK(edge_count == new_edge_count);
        //next: check degrees stay the same.
        OK (LAGraph_Cached_OutDegree (G_new, msg)) ;

        OK (LAGraph_Vector_IsEqual (
            &ok, G->out_degree, G_new->out_degree, msg)) ;
        TEST_CHECK (ok) ;

        OK (LAGraph_Delete (&G, msg)) ;
        OK (LAGraph_Delete (&G_new, msg)) ;
    }

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------
    LAGraph_Finalize (msg) ;
    #endif
}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"SwapEdges", test_SwapEdges},   
    {NULL, NULL}
} ;
