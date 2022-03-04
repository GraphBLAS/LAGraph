//------------------------------------------------------------------------------
// LAGraph/src/test/test_Property_Degree.c:  test LAGraph_Property_*Degree
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
char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL ;
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
// check_degree: check a row or column degree vector
//------------------------------------------------------------------------------

void check_degree
(
    GrB_Vector Degree,
    GrB_Index n,
    const int *degree
)
{
    GrB_Index n2 ;
    OK (GrB_Vector_size (&n2, Degree)) ;
    TEST_CHECK (n == n2) ;
    for (int k = 0 ; k < n ; k++)
    {
        int64_t degk ;
        GrB_Info info = GrB_Vector_extractElement (&degk, Degree, k) ;
        TEST_CHECK (info == GrB_NO_VALUE || info == GrB_SUCCESS) ;
        if (info == GrB_NO_VALUE)
        {
            TEST_CHECK (degree [k] == 0) ;
        }
        else
        {
            TEST_CHECK (degree [k] == degk) ;
        }
    }
}

//------------------------------------------------------------------------------
// test_Property_Degree:  test LAGraph_Property_*Degree
//------------------------------------------------------------------------------

typedef struct
{
    const char *name ;
    const int rowdeg [67] ;
    const int coldeg [67] ;
}
matrix_info ;

const matrix_info files [ ] =
{
    { "A.mtx", 
        { 3, 5, 5, 5, 3, 4, 5,  },
        { 3, 5, 5, 5, 3, 4, 5,  }, },
     { "LFAT5.mtx", 
        { 3, 2, 2, 4, 4, 3, 3, 5, 5, 2, 2, 4, 4, 3,  },
        { 3, 2, 2, 4, 4, 3, 3, 5, 5, 2, 2, 4, 4, 3,  }, },
     { "cover.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "cover_structure.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "full.mtx", 
        { 3, 3, 3,  },
        { 3, 3, 3,  }, },
     { "full_symmetric.mtx", 
        { 4, 4, 4, 4,  },
        { 4, 4, 4, 4,  }, },
     { "karate.mtx", 
        { 16, 9, 10, 6, 3, 4, 4, 4, 5, 2, 3, 1, 2, 5, 2, 2, 2, 2, 2, 3,
          2, 2, 2, 5, 3, 3, 2, 4, 3, 4, 4, 6, 12, 17,  },
        { 16, 9, 10, 6, 3, 4, 4, 4, 5, 2, 3, 1, 2, 5, 2, 2, 2, 2, 2, 3,
          2, 2, 2, 5, 3, 3, 2, 4, 3, 4, 4, 6, 12, 17,  }, },
     { "ldbc-cdlp-directed-example.mtx", 
        { 3, 2, 2, 2, 3, 2, 3, 1,  },
        { 2, 2, 2, 1, 3, 4, 3, 1,  }, },
     { "ldbc-cdlp-undirected-example.mtx", 
        { 3, 2, 2, 3, 4, 3, 3, 4,  },
        { 3, 2, 2, 3, 4, 3, 3, 4,  }, },
     { "ldbc-directed-example-bool.mtx", 
        { 2, 3, 4, 0, 3, 2, 1, 1, 1, 0,  },
        { 2, 0, 3, 5, 3, 0, 0, 2, 0, 2,  }, },
     { "ldbc-directed-example-unweighted.mtx", 
        { 2, 3, 4, 0, 3, 2, 1, 1, 1, 0,  },
        { 2, 0, 3, 5, 3, 0, 0, 2, 0, 2,  }, },
     { "ldbc-directed-example.mtx", 
        { 2, 3, 4, 0, 3, 2, 1, 1, 1, 0,  },
        { 2, 0, 3, 5, 3, 0, 0, 2, 0, 2,  }, },
     { "ldbc-undirected-example-bool.mtx", 
        { 2, 4, 2, 3, 5, 2, 3, 2, 1,  },
        { 2, 4, 2, 3, 5, 2, 3, 2, 1,  }, },
     { "ldbc-undirected-example-unweighted.mtx", 
        { 2, 4, 2, 3, 5, 2, 3, 2, 1,  },
        { 2, 4, 2, 3, 5, 2, 3, 2, 1,  }, },
     { "ldbc-undirected-example.mtx", 
        { 2, 4, 2, 3, 5, 2, 3, 2, 1,  },
        { 2, 4, 2, 3, 5, 2, 3, 2, 1,  }, },
     { "ldbc-wcc-example.mtx", 
        { 3, 3, 5, 5, 5, 2, 1, 3, 1, 2,  },
        { 3, 3, 5, 5, 5, 2, 1, 3, 1, 2,  }, },
     { "matrix_bool.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_fp32.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_fp32_structure.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_fp64.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_int16.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_int32.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_int64.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_int8.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_uint16.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_uint32.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_uint64.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "matrix_uint8.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "msf1.mtx", 
        { 2, 2, 1, 1, 1, 1,  },
        { 1, 1, 2, 2, 0, 2,  }, },
     { "msf2.mtx", 
        { 2, 3, 3, 2, 1, 1, 0, 0,  },
        { 0, 1, 1, 1, 2, 2, 2, 3,  }, },
     { "msf3.mtx", 
        { 2, 2, 2, 1, 0,  },
        { 0, 1, 1, 2, 3,  }, },
     { "structure.mtx", 
        { 2, 2, 1, 2, 1, 1, 3,  },
        { 1, 1, 3, 2, 2, 2, 1,  }, },
     { "sample.mtx", 
        { 3, 2, 1, 2, 2, 1, 1, 0,  },
        { 0, 1, 3, 1, 3, 1, 1, 2,  }, },
     { "sample2.mtx", 
        { 2, 3, 4, 3, 5, 5, 3, 3,  },
        { 2, 3, 4, 3, 5, 5, 3, 3,  }, },
     { "skew_fp32.mtx", 
        { 3, 3, 3, 4, 3, 4,  },
        { 3, 3, 3, 4, 3, 4,  }, },
     { "skew_fp64.mtx", 
        { 3, 3, 3, 4, 3, 4,  },
        { 3, 3, 3, 4, 3, 4,  }, },
     { "skew_int16.mtx", 
        { 3, 3, 3, 4, 3, 4,  },
        { 3, 3, 3, 4, 3, 4,  }, },
     { "skew_int32.mtx", 
        { 3, 3, 3, 4, 3, 4,  },
        { 3, 3, 3, 4, 3, 4,  }, },
     { "skew_int64.mtx", 
        { 3, 3, 3, 4, 3, 4,  },
        { 3, 3, 3, 4, 3, 4,  }, },
     { "skew_int8.mtx", 
        { 3, 3, 3, 4, 3, 4,  },
        { 3, 3, 3, 4, 3, 4,  }, },
     { "tree-example.mtx", 
        { 1, 1, 2, 3, 2, 1,  },
        { 1, 1, 2, 3, 2, 1,  }, },
     { "west0067.mtx", 
        { 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5,
          3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 3, 3, 3, 3, 5,
          5, 5, 5, 5, 6, 3, 3, 3, 3, 4, 4, 4, 4, 4, 6, 1, 5, 5, 5, 5,
          5, 5, 5, 5, 5, 5, 5,  },
        { 10, 4, 4, 4, 4, 3, 5, 3, 3, 3, 3, 2, 5, 5, 5, 5, 4, 5, 2, 10,
          3, 3, 3, 3, 3, 4, 4, 4, 4, 3, 10, 3, 3, 3, 3, 3, 10, 5, 5, 5,
          5, 4, 5, 4, 4, 4, 4, 3, 10, 3, 3, 3, 3, 3, 10, 5, 5, 5, 5, 4,
          5, 4, 4, 4, 4, 3, 5,  }, },
     { "west0067_jumbled.mtx", 
        { 3, 3, 3, 3, 5, 5, 5, 5, 5, 6, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5,
          3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 3, 3, 3, 3, 5,
          5, 5, 5, 5, 6, 3, 3, 3, 3, 4, 4, 4, 4, 4, 6, 1, 5, 5, 5, 5,
          5, 5, 5, 5, 5, 5, 5,  },
        { 10, 4, 4, 4, 4, 3, 5, 3, 3, 3, 3, 2, 5, 5, 5, 5, 4, 5, 2, 10,
          3, 3, 3, 3, 3, 4, 4, 4, 4, 3, 10, 3, 3, 3, 3, 3, 10, 5, 5, 5,
          5, 4, 5, 4, 4, 4, 4, 3, 10, 3, 3, 3, 3, 3, 10, 5, 5, 5, 5, 4,
          5, 4, 4, 4, 4, 3, 5,  }, },
    { "", { 0 }, { 0 }}
} ;

//-----------------------------------------------------------------------------
// test_Property_Degree
//-----------------------------------------------------------------------------

void test_Property_Degree (void)
{
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        const int *rowdeg = files [k].rowdeg ;
        const int *coldeg = files [k].coldeg ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // construct the graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        for (int trial = 0 ; trial <= 2 ; trial++)
        {
            // create the G->rowdegree property and check it
            OK (LAGraph_Property_RowDegree (G, msg)) ;
            GrB_Index n ;
            OK (GrB_Matrix_nrows (&n, G->A)) ;
            check_degree (G->rowdegree, n, rowdeg) ;

            if (trial == 2)
            {
                // use G->AT to compute G->coldegree 
                OK (LAGraph_DeleteProperties (G, msg)) ;
                OK (LAGraph_Property_AT (G, msg)) ;
            }

            // create the G->ColDegree property and check it
            OK (LAGraph_Property_ColDegree (G, msg)) ;
            OK (GrB_Matrix_ncols (&n, G->A)) ;
            check_degree (G->coldegree, n, coldeg) ;
        }

        OK (LAGraph_Delete (&G, msg)) ;
    }

    // check error handling
    int status = LAGraph_Property_RowDegree (NULL, msg) ;
    printf ("\nstatus: %d, msg: %s\n", status, msg) ;
    TEST_CHECK (status == GrB_NULL_POINTER) ;
    status = LAGraph_Property_ColDegree (NULL, msg) ;
    printf ("status: %d, msg: %s\n", status, msg) ;
    TEST_CHECK (status == GrB_NULL_POINTER) ;

    teardown ( ) ;
}

//-----------------------------------------------------------------------------
// test_Property_Degree_brutal
//-----------------------------------------------------------------------------

#if LAGRAPH_SUITESPARSE
void test_Property_Degree_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        const int *rowdeg = files [k].rowdeg ;
        const int *coldeg = files [k].coldeg ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // construct the graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        for (int trial = 0 ; trial <= 2 ; trial++)
        {
            // create the G->rowdegree property and check it
            LG_BRUTAL (LAGraph_Property_RowDegree (G, msg)) ;
            GrB_Index n ;
            OK (GrB_Matrix_nrows (&n, G->A)) ;
            check_degree (G->rowdegree, n, rowdeg) ;

            if (trial == 2)
            {
                // use G->AT to compute G->coldegree 
                OK (LAGraph_DeleteProperties (G, msg)) ;
                OK (LAGraph_Property_AT (G, msg)) ;
            }

            // create the G->ColDegree property and check it
            LG_BRUTAL (LAGraph_Property_ColDegree (G, msg)) ;
            OK (GrB_Matrix_ncols (&n, G->A)) ;
            check_degree (G->coldegree, n, coldeg) ;
        }

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
    { "Property_Degree", test_Property_Degree },
    #if LAGRAPH_SUITESPARSE
    { "Property_Degree_brutal", test_Property_Degree_brutal },
    #endif
    { NULL, NULL }
} ;

