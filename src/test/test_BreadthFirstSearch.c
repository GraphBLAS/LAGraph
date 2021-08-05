//----------------------------------------------------------------------------
// LAGraph/src/test/test_BreadthFirstSearch.cpp: test cases for triangle
// counting algorithms
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
#include <graph_zachary_karate.h>

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

GrB_Index const SRC = 30;
GrB_Index const LEVELS30[] = {2, 1, 2, 2, 3, 3, 3, 2, 1, 2, 3, 3,
                              3, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2,
                              3, 3, 2, 2, 2, 2, 0, 2, 1, 1};
GrB_Index const PARENT30[] = { 1, 30,  1,  1,  0,  0,  0,  1, 30, 33,  0,  0,
                               0,  1, 32, 32,  5,  1, 32,  1, 32,  1, 32, 32,
                              27, 23, 33, 33, 33, 32, 30, 32, 30, 30};

#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    LAGraph_Kind kind ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    { LAGRAPH_ADJACENCY_UNDIRECTED, "A.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "cover.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "jagmesh7.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "ldbc-cdlp-directed-example.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-cdlp-undirected-example.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "ldbc-directed-example.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-undirected-example.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-wcc-example.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "LFAT5.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "msf1.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "msf2.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "msf3.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "sample2.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "sample.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "olm1000.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "bcsstk13.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "cryg2500.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "tree-example.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "west0067.mtx" }, 
    { LAGRAPH_ADJACENCY_UNDIRECTED, "karate.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_bool.mtx" }, 
    { LAGRAPH_ADJACENCY_DIRECTED,   "skew_fp32.mtx" }, 
    { LAGRAPH_UNKNOWN, "" }, 
} ;

//****************************************************************************
bool check_karate_parents30(GrB_Vector parents)
{
    // TODO: this may not work in multithreaded code (w/ benign races)
    GrB_Index n = 0;
    TEST_CHECK(0 == GrB_Vector_size(&n, parents));
    TEST_CHECK(ZACHARY_NUM_NODES == n);
    TEST_CHECK(0 == GrB_Vector_nvals(&n, parents));
    TEST_CHECK(ZACHARY_NUM_NODES == n);

    int64_t parent_id;
    for (GrB_Index ix = 0; ix < ZACHARY_NUM_NODES; ++ix)
    {
        TEST_CHECK(0 == GrB_Vector_extractElement(&parent_id, parents, ix));
        TEST_CHECK(parent_id == PARENT30[ix]);
        TEST_MSG("Parent check failed for node %ld: ans,comp = %ld,%ld\n",
                 ix, PARENT30[ix], parent_id);
    }

    return true;
}

//****************************************************************************
bool check_karate_levels30(GrB_Vector levels)
{
    GrB_Index n = 0;
    TEST_CHECK(0 == GrB_Vector_size(&n, levels) );
    TEST_CHECK(ZACHARY_NUM_NODES == n);
    TEST_CHECK(0 == GrB_Vector_nvals(&n, levels) );
    TEST_CHECK(ZACHARY_NUM_NODES == n);

    int64_t lvl;
    for (GrB_Index ix = 0; ix < ZACHARY_NUM_NODES; ++ix)
    {
        TEST_CHECK(0 == GrB_Vector_extractElement(&lvl, levels, ix) );
        TEST_CHECK(lvl == LEVELS30[ix] );
        TEST_MSG("Level check failed for node %ld: ans,comp = %ld,%ld\n",
                 ix, LEVELS30[ix], lvl);
    }

    return true;
}

//****************************************************************************
void setup(void)
{
    LAGraph_Init(msg);
    int retval;
    GrB_Matrix A = NULL;

    TEST_CHECK(0 == GrB_Matrix_new(&A, GrB_UINT32,
                                   ZACHARY_NUM_NODES, ZACHARY_NUM_NODES) );
    TEST_CHECK(0 == GrB_Matrix_build(A, ZACHARY_I, ZACHARY_J, ZACHARY_V,
                                     ZACHARY_NUM_EDGES, GrB_LOR) );

    retval = LAGraph_New(&G, &A, GrB_UINT32, LAGRAPH_ADJACENCY_UNDIRECTED, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
}

//****************************************************************************
void teardown(void)
{
    int retval = LAGraph_Delete(&G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    G = NULL;
    LAGraph_Finalize(msg);
}

//****************************************************************************
void test_BreadthFirstSearch_invalid_graph(void)
{
    setup();
    int retval;
    LAGraph_Graph graph = NULL;

    retval = LAGraph_BreadthFirstSearch(NULL, NULL, graph, 0, false, msg);
    TEST_CHECK(retval == -101);
    TEST_MSG("retval = %d (%s)", retval, msg);

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_invalid_src(void)
{
    setup();
    int retval;
    GrB_Index n;
    TEST_CHECK(0 == GrB_Matrix_nrows(&n, (G->A)));

    GrB_Vector parent    = NULL;

    retval = LAGraph_BreadthFirstSearch(NULL, NULL, G, n, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch(NULL, &parent, G, n, false, msg);
    TEST_CHECK(retval == -102);
    TEST_MSG("retval = %d (%s)", retval, msg);

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_neither(void)
{
    setup();
    int retval;

    // Suitesparse extensions
    retval = LAGraph_BreadthFirstSearch(NULL, NULL, G, 0, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    // Suitesparse extensions
    retval = LAGraph_BreadthFirstSearch(NULL, NULL, G, 0, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_parent(void)
{
    setup();
    int retval;

    GrB_Vector parent    = NULL;
    GrB_Vector parent_do = NULL;

    retval = LAGraph_BreadthFirstSearch(NULL, &parent, G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_parents30(parent));

    retval = LG_check_bfs (NULL, parent, G, 30, msg) ;
    TEST_CHECK (retval == 0) ;

    //printf("parent BFS (!pushpull) source = 30\n");
    //LAGraph_Vector_print_type(parent, GrB_INT32, 3, stdout, msg);

    retval = LAGraph_BreadthFirstSearch(NULL, &parent_do, G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_parents30(parent_do));

    retval = LG_check_bfs (NULL, parent_do, G, 30, msg) ;
    TEST_CHECK (retval == 0) ;

    //printf("parent BFS (pushpull) source = 30\n");
    //LAGraph_Vector_print_type(parent_do, GrB_INT32, 3, stdout, msg);

    TEST_CHECK(0 == GrB_free(&parent));
    TEST_CHECK(0 == GrB_free(&parent_do));

    GrB_Index n = 0 ;
    TEST_CHECK (0 == GrB_Matrix_nrows (&n, G->A)) ;
    for (int pushpull = 0 ; pushpull <= 1 ; pushpull++)
    {
        for (GrB_Index src = 0 ; src < n ; src++)
        {
            retval = LAGraph_BreadthFirstSearch(NULL, &parent,
                                                G, src, (bool) pushpull, msg);
            TEST_CHECK(retval == 0);
            retval = LG_check_bfs (NULL, parent, G, src, msg) ;
            TEST_CHECK (retval == 0) ;
            TEST_CHECK(0 == GrB_free(&parent));
        }
    }

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_level(void)
{
    setup();
    int retval;

    GrB_Vector level    = NULL;
    GrB_Vector level_do = NULL;

    retval = LAGraph_BreadthFirstSearch(&level, NULL, G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_levels30(level));

    retval = LG_check_bfs (level, NULL, G, 30, msg) ;
    TEST_CHECK (retval == 0) ;

    TEST_CHECK(0 == GrB_free(&level));

    retval = LAGraph_BreadthFirstSearch(&level_do, NULL, G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_levels30(level_do));

    retval = LG_check_bfs (level_do, NULL, G, 30, msg) ;
    TEST_CHECK (retval == 0) ;

    TEST_CHECK(0 == GrB_free(&level_do));

    GrB_Index n = 0 ;
    TEST_CHECK (0 == GrB_Matrix_nrows (&n, G->A)) ;
    for (int pushpull = 0 ; pushpull <= 1 ; pushpull++)
    {
        for (GrB_Index src = 0 ; src < n ; src++)
        {
            retval = LAGraph_BreadthFirstSearch(&level, NULL,
                                                G, src, (bool) pushpull, msg);
            TEST_CHECK(retval == 0);
            retval = LG_check_bfs (level, NULL, G, src, msg) ;
            TEST_CHECK (retval == 0) ;
            TEST_CHECK(0 == GrB_free(&level));
        }
    }

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_both(void)
{
    setup();
    int retval;

    GrB_Vector parent    = NULL;
    GrB_Vector level    = NULL;
    retval = LAGraph_BreadthFirstSearch(&level, &parent,
                                        G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_levels30(level));
    TEST_CHECK(check_karate_parents30(parent));

    retval = LG_check_bfs (level, parent, G, 30, msg) ;
    TEST_CHECK (retval == 0) ;

    TEST_CHECK(0 == GrB_free(&parent));
    TEST_CHECK(0 == GrB_free(&level));

    GrB_Vector parent_do = NULL;
    GrB_Vector level_do = NULL;
    retval = LAGraph_BreadthFirstSearch(&level_do, &parent_do,
                                        G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_levels30(level_do));
    TEST_CHECK(check_karate_parents30(parent_do));
    retval = LG_check_bfs (level_do, parent_do, G, 30, msg) ;
    TEST_CHECK (retval == 0) ;
    TEST_CHECK(0 == GrB_free(&parent_do));
    TEST_CHECK(0 == GrB_free(&level_do));

    GrB_Index n = 0 ;
    TEST_CHECK (0 == GrB_Matrix_nrows (&n, G->A)) ;
    for (int pushpull = 0 ; pushpull <= 1 ; pushpull++)
    {
        for (GrB_Index src = 0 ; src < n ; src++)
        {
            retval = LAGraph_BreadthFirstSearch(&level, &parent,
                                                G, src, (bool) pushpull, msg);
            TEST_CHECK(retval == 0);
            retval = LG_check_bfs (level, parent, G, src, msg) ;
            TEST_CHECK (retval == 0) ;
            TEST_CHECK(0 == GrB_free(&parent));
            TEST_CHECK(0 == GrB_free(&level));
        }
    }

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_brutal(void)
{
    LAGraph_Init(msg);
    GrB_Matrix A = NULL ;
    GrB_Type atype = NULL ;

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        LAGraph_Kind kind = files [k].kind ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\nMatrix: %s\n", aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // create its properties
        OK (LAGraph_Property_AT (G, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;

        OK (LAGraph_Property_RowDegree (G, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;

        OK (LAGraph_Property_ColDegree (G, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;

        GrB_Index n = 0 ;
        OK (GrB_Matrix_nrows (&n, G->A)) ;

        // run the BFS
        int64_t step = (n > 100) ? (3*n/4) : ((n/4) + 1) ;
        for (int64_t src = 0 ; src < n ; src += step)
        {
            printf ("src %ld\n", src) ;
            GrB_Vector parent = NULL ;
            GrB_Vector level = NULL ;
            for (int pushpull = 0 ; pushpull <= 1 ; pushpull++)
            {

                OK (LAGraph_BreadthFirstSearch (&level, &parent,
                    G, src, (bool) pushpull, msg)) ;
                OK (LG_check_bfs (level, parent, G, src, msg)) ;
                int64_t maxlevel, nvisited ;
                OK (GrB_reduce (&maxlevel, NULL, GrB_MAX_MONOID_INT64,
                    level, NULL)) ;
                OK (GrB_Vector_nvals (&nvisited, level)) ;
                printf (" n: %ld max level: %ld nvisited: %ld\n", n,
                    maxlevel, nvisited) ;

                OK (GrB_free(&parent));
                OK (GrB_free(&level));

                OK (LAGraph_BreadthFirstSearch (NULL, &parent,
                    G, src, (bool) pushpull, msg)) ;
                OK (LG_check_bfs (NULL, parent, G, src, msg)) ;
                OK (GrB_free(&parent));

                OK (LAGraph_BreadthFirstSearch (&level, NULL,
                    G, src, (bool) pushpull, msg)) ;
                OK (LG_check_bfs (level, NULL, G, src, msg)) ;
                OK (GrB_free(&level));
            }
        }
    }

    LAGraph_Finalize(msg);
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"BreadthFirstSearch_invalid_graph", test_BreadthFirstSearch_invalid_graph},
    {"BreadthFirstSearch_invalid_src", test_BreadthFirstSearch_invalid_src},
    {"BreadthFirstSearch_neither", test_BreadthFirstSearch_neither},
    {"BreadthFirstSearch_parent", test_BreadthFirstSearch_parent},
    {"BreadthFirstSearch_level", test_BreadthFirstSearch_level},
    {"BreadthFirstSearch_both", test_BreadthFirstSearch_both},
    {"BreadthFirstSearch_brutal", test_BreadthFirstSearch_brutal},
    {NULL, NULL}
};

