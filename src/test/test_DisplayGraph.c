//------------------------------------------------------------------------------
// LAGraph/src/test/test_DisplayGraph.c:  test LAGraph_DisplayGraph
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
// prwhat: print what should be printed
//------------------------------------------------------------------------------

const char *prwhat (int pr)
{
    switch (pr)
    {
        case -1: return ("nothing") ;
        case  0: return ("single line") ;
        case  1: return ("terse") ;
        case  2: return ("summary") ;
        case  3: return ("all") ;
        case  4: return ("summary (doubles in full precision)") ;
        case  5: return ("all (doubles in full precision)") ;
        default: ;
    }
    return (NULL) ;
}

//------------------------------------------------------------------------------
// test_DisplayGraph:  test LAGraph_DisplayGraph
//------------------------------------------------------------------------------

typedef struct
{
    LAGraph_Kind kind ;
    int ndiag ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    LAGRAPH_ADJACENCY_DIRECTED,   0, "cover.mtx",
    LAGRAPH_ADJACENCY_DIRECTED,   0, "ldbc-directed-example.mtx",
    LAGRAPH_ADJACENCY_UNDIRECTED, 0, "ldbc-undirected-example.mtx",
    LAGRAPH_ADJACENCY_DIRECTED,   2, "west0067.mtx",
    LAGRAPH_UNKNOWN,              0, ""
} ;

void test_DisplayGraph (void)
{
    setup ( ) ;

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        LAGraph_Kind kind = files [k].kind ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        if (atype == GrB_FP64)
        {
            OK (GrB_Matrix_setElement (A, 3.14159265358979323, 0, 1)) ;
        }

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // display the graph
        for (int trial = 0 ; trial <= 1 ; trial++)
        {
            printf ("\n############################# TRIAL: %d\n", trial) ;
            for (int pr = -1 ; pr <= 5 ; pr++)
            {
                printf ("\n########### %s: pr: %d (%s)\n",
                    aname, pr, prwhat (pr)) ;
                OK (LAGraph_DisplayGraph (G, pr, stdout, msg)) ;
            }
            OK (LAGraph_Property_AT (G, msg)) ;
            OK (LAGraph_Property_ASymmetricStructure (G, msg)) ;
            OK (LAGraph_Property_NDiag (G, msg)) ;
            TEST_CHECK (G->ndiag == files [k].ndiag) ;
        }

        // free the graph
        OK (LAGraph_Delete (&G, msg)) ;
        TEST_CHECK (G == NULL) ;
    }

    TEST_CHECK (prwhat (999) == NULL) ;
    teardown ( ) ;
}

//------------------------------------------------------------------------------
// test_DisplayGraph_failures:  test error handling of LAGraph_DisplayGraph
//------------------------------------------------------------------------------

#if 0
void test_DisplayGraph_failures (void)      // TODO
{
    setup ( ) ;

    // G cannot be NULL
    TEST_CHECK (LAGraph_New (NULL, NULL, NULL, 0, msg) == -1) ;
    printf ("\nmsg: %s\n", msg) ;

    // create a graph with no adjacency matrix; this is OK, since the intent is
    // to create a graph for which the adjacency matrix can be defined later,
    // via assigning it to G->A.  However, the graph will be declared invalid
    // by LAGraph_CheckGraph since G->A is NULL.
    OK (LAGraph_New (&G, NULL, NULL, 0, msg)) ;
    TEST_CHECK (LAGraph_CheckGraph (G, msg) == -2) ;
    printf ("msg: %s\n", msg) ;
    OK (LAGraph_Delete (&G, msg)) ;
    TEST_CHECK (G == NULL) ;
    OK (LAGraph_Delete (&G, msg)) ;
    TEST_CHECK (G == NULL) ;
    OK (LAGraph_Delete (NULL, msg)) ;
    teardown ( ) ;
}
#endif

//-----------------------------------------------------------------------------
// test_DisplayGraph_brutal
//-----------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_DisplayGraph_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        LAGraph_Kind kind = files [k].kind ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        if (atype == GrB_FP64)
        {
            OK (GrB_Matrix_setElement (A, 3.14159265358979323, 0, 1)) ;
        }
        OK (GrB_wait (A, GrB_MATERIALIZE)) ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;

        // display the graph
        for (int trial = 0 ; trial <= 1 ; trial++)
        {
            printf ("\n############################# TRIAL: %d\n", trial) ;
            for (int pr = -1 ; pr <= 5 ; pr++)
            {
                printf ("\n########### %s: pr: %d (%s)\n",
                    aname, pr, prwhat (pr)) ;
                if (pr == 3 || pr == 5)
                {
                    printf ("skipped for brutal tests\n") ;
                }
                else
                {
                    LG_BRUTAL (LAGraph_DisplayGraph (G, pr, stdout, msg)) ;
                }
            }
            OK (LAGraph_Property_AT (G, msg)) ;
            OK (LAGraph_Property_ASymmetricStructure (G, msg)) ;
            OK (LAGraph_Property_NDiag (G, msg)) ;
        }

        // free the graph
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
    { "DisplayGraph", test_DisplayGraph },
//  { "DisplayGraph_failures", test_DisplayGraph_failures },
    { "DisplayGraph_brutal", test_DisplayGraph_brutal },
    { NULL, NULL }
} ;

