//----------------------------------------------------------------------------
// LAGraph/experimental/test/test_louvain.c: test cases for LAGraph_louvain
//----------------------------------------------------------------------------

// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Roi Lipman and Gabriel Gomez, FalkorDB

//-----------------------------------------------------------------------------

// Smoke test for the Louvain community-detection method: run the algorithm on
// a few small graphs and confirm the result is a well-formed community vector.

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix    A   = NULL ;
LAGraph_Graph G   = NULL ;
GrB_Vector    com = NULL ;

#define LEN 512
char filename [LEN+1] ;

// Louvain parameters used throughout the smoke test.
#define ITERMAX  6      // max modularity-improvement sweeps per level
#define LEVELMAX 2      // max improve-and-condense levels
#define EPSILON  1e-5f  // min modularity change considered an improvement

typedef struct
{
    const char *name ;
}
matrix_info ;

// pattern / boolean symmetric matrices: LAGraph_louvain requires G->A to be a
// square boolean adjacency matrix, which is exactly what these files load as.
const matrix_info files [ ] =
{
    { "A.mtx" },        // tiny 7x7 graph
    { "karate.mtx" },   // Zachary's karate club, the classic clustering example
    { "" },
} ;

//****************************************************************************
void test_louvain (void)
{
    OK (LAGraph_Init (msg)) ;

    for (int k = 0 ; ; k++)
    {
        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0)
		{
			break ;
		}

        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;

        // construct an undirected graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // Louvain operates on a simple graph
        OK (LAGraph_DeleteSelfEdges (G, msg)) ;

        GrB_Index n ;
        OK (GrB_Matrix_nrows (&n, G->A)) ;

        // compute the communities with LAGraph_louvain
        int result = LAGraph_louvain (&com, G, ITERMAX, LEVELMAX, EPSILON, msg) ;
        TEST_CHECK (result == GrB_SUCCESS) ;

        //----------------------------------------------------------------------
        // check the result is a well-formed community vector
        //----------------------------------------------------------------------

        // one community id per node
        TEST_CHECK (com != NULL) ;
        GrB_Index com_size ;
        OK (GrB_Vector_size (&com_size, com)) ;
        TEST_CHECK (com_size == n) ;
        TEST_MSG ("community vector size %g, expected %g",
            (double) com_size, (double) n) ;

        // the result must be a full vector: one community id per node
        GrB_Index nvals ;
        OK (GrB_Vector_nvals (&nvals, com)) ;
        TEST_CHECK (nvals == n) ;
        TEST_MSG ("community vector has %g entries, expected %g",
            (double) nvals, (double) n) ;

        // print the result
        LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT ;
        printf ("\nlouvain (computed communities):\n") ;
        OK (LAGraph_Vector_Print (com, pr, stdout, msg)) ;

        // community ids are node indices, so the largest id must be < n;
        // ids are unsigned, so this bounds every id to the valid range [0, n)
        uint64_t max_id = 0 ;
        OK (GrB_reduce (&max_id, NULL, GrB_MAX_MONOID_UINT64, com, NULL)) ;
        TEST_CHECK (max_id < n) ;
        TEST_MSG ("largest community id %g, must be < %g",
            (double) max_id, (double) n) ;

        OK (GrB_free (&com)) ;
        OK (LAGraph_Delete (&G, msg)) ;
    }

    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_louvain_errors (void)
{
    OK (LAGraph_Init (msg)) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    TEST_MSG ("Loading of adjacency matrix failed") ;
    OK (fclose (f)) ;

    // construct an undirected graph G with adjacency matrix A
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;
    OK (LAGraph_DeleteSelfEdges (G, msg)) ;

    // com is NULL
    int result = LAGraph_louvain (NULL, G, ITERMAX, LEVELMAX, EPSILON, msg) ;
    printf ("\nresult: %d\n", result) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // G is NULL
    result = LAGraph_louvain (&com, NULL, ITERMAX, LEVELMAX, EPSILON, msg) ;
    printf ("result: %d\n", result) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
}

//****************************************************************************

TEST_LIST = {
    {"louvain", test_louvain},
    {"louvain_errors", test_louvain_errors},
    {NULL, NULL}
} ;
