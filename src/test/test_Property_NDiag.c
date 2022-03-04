//------------------------------------------------------------------------------
// LAGraph/src/test/test_Property_NDiag.c:  test LAGraph_Property_NDiag
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LAGraph_test.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

LAGraph_Graph G = NULL ;
GrB_Info info ;
char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL ;
#define LEN 512
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// test matrices
//------------------------------------------------------------------------------

typedef struct
{
    GrB_Index ndiag ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = 
{
     0, "A.mtx",
    14, "LFAT5.mtx",
  2003, "bcsstk13.mtx",
//   3, "complex.mtx",
     0, "cover.mtx",
     0, "cover_structure.mtx",
  2500, "cryg2500.mtx",
     3, "full.mtx",
     4, "full_symmetric.mtx",
  1138, "jagmesh7.mtx",
     0, "karate.mtx",
     0, "ldbc-cdlp-directed-example.mtx",
     0, "ldbc-cdlp-undirected-example.mtx",
     0, "ldbc-directed-example-bool.mtx",
     0, "ldbc-directed-example-unweighted.mtx",
     0, "ldbc-directed-example.mtx",
     0, "ldbc-undirected-example-bool.mtx",
     0, "ldbc-undirected-example-unweighted.mtx",
     0, "ldbc-undirected-example.mtx",
     0, "ldbc-wcc-example.mtx",
     0, "matrix_bool.mtx",
     0, "matrix_fp32.mtx",
     0, "matrix_fp32_structure.mtx",
     0, "matrix_fp64.mtx",
     0, "matrix_int16.mtx",
     0, "matrix_int32.mtx",
     0, "matrix_int64.mtx",
     0, "matrix_int8.mtx",
     0, "matrix_uint16.mtx",
     0, "matrix_uint32.mtx",
     0, "matrix_uint64.mtx",
     0, "matrix_uint8.mtx",
     0, "msf1.mtx",
     0, "msf2.mtx",
     0, "msf3.mtx",
  1000, "olm1000.mtx",
     0, "structure.mtx",
     0, "sample.mtx",
     0, "sample2.mtx",
     0, "skew_fp32.mtx",
     0, "skew_fp64.mtx",
     0, "skew_int16.mtx",
     0, "skew_int32.mtx",
     0, "skew_int64.mtx",
     0, "skew_int8.mtx",
     0, "tree-example.mtx",
     2, "west0067.mtx",
     2, "west0067_jumbled.mtx",
     0, ""
} ;

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
// test_Property_NDiag:  test LAGraph_Property_NDiag
//------------------------------------------------------------------------------

void test_Property_NDiag (void)
{

    //--------------------------------------------------------------------------
    // start up the test
    //--------------------------------------------------------------------------

    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth file
        //----------------------------------------------------------------------

        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;

        //----------------------------------------------------------------------
        // construct a directed graph and count self-edges
        //----------------------------------------------------------------------

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        OK (LAGraph_Property_NDiag (G, msg)) ;
        TEST_CHECK (G->ndiag == files [k].ndiag) ;

        //----------------------------------------------------------------------
        // free the graph
        //----------------------------------------------------------------------

        OK (LAGraph_Delete (&G, msg)) ;
    }

    TEST_CHECK (LAGraph_Property_NDiag (NULL, msg) == GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // finish the test
    //--------------------------------------------------------------------------

    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_Property_NDiag_brutal
//------------------------------------------------------------------------------

#if LAGRAPH_SUITESPARSE
void test_Property_NDiag_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    for (int k = 0 ; ; k++)
    {

        //----------------------------------------------------------------------
        // load in the kth file
        //----------------------------------------------------------------------

        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Failed to load %s\n", aname) ;

        //----------------------------------------------------------------------
        // construct a directed graph and count self-edges
        //----------------------------------------------------------------------

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        LG_BRUTAL (LAGraph_Property_NDiag (G, msg)) ;
        TEST_CHECK (G->ndiag == files [k].ndiag) ;

        //----------------------------------------------------------------------
        // free the graph
        //----------------------------------------------------------------------

        OK (LAGraph_Delete (&G, msg)) ;
    }

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "NDiag", test_Property_NDiag },
    #if LAGRAPH_SUITESPARSE
    { "NDiag_brutal", test_Property_NDiag_brutal },
    #endif
    { NULL, NULL }
} ;

