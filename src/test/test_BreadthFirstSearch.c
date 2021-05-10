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

    retval = LAGraph_BreadthFirstSearch_vanilla(NULL, NULL, graph, 0, false, msg);
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

    retval = LAGraph_BreadthFirstSearch_vanilla(NULL, NULL, G, n, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch_vanilla(NULL, &parent, G, n, false, msg);
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

    // no extensions
    retval = LAGraph_BreadthFirstSearch_vanilla(NULL, NULL, G, 0, true, msg);
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
    GrB_Vector parent_v  = NULL;

    retval = LAGraph_BreadthFirstSearch(NULL, &parent, G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch(NULL, &parent_do, G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch_vanilla(NULL, &parent_v, G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK(0 == GrB_free(&parent));
    TEST_CHECK(0 == GrB_free(&parent_do));
    TEST_CHECK(0 == GrB_free(&parent_v));

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_level(void)
{
    setup();
    int retval;

    GrB_Vector level    = NULL;
    GrB_Vector level_do = NULL;
    GrB_Vector level_v  = NULL;

    retval = LAGraph_BreadthFirstSearch(&level, NULL, G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch(&level_do, NULL, G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch_vanilla(&level_v, NULL, G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK(0 == GrB_free(&level));
    TEST_CHECK(0 == GrB_free(&level_do));
    TEST_CHECK(0 == GrB_free(&level_v));

    teardown();
}

//****************************************************************************
void test_BreadthFirstSearch_both(void)
{
    setup();
    int retval;

    GrB_Vector parent    = NULL;
    GrB_Vector parent_do = NULL;
    GrB_Vector parent_v  = NULL;
    GrB_Vector level    = NULL;
    GrB_Vector level_do = NULL;
    GrB_Vector level_v  = NULL;

    retval = LAGraph_BreadthFirstSearch(&level, &parent,
                                        G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch(&level_do, &parent_do,
                                        G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_BreadthFirstSearch_vanilla(&level_v, &parent_v,
                                                G, 30, false, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK(0 == GrB_free(&parent));
    TEST_CHECK(0 == GrB_free(&parent_do));
    TEST_CHECK(0 == GrB_free(&parent_v));
    TEST_CHECK(0 == GrB_free(&level));
    TEST_CHECK(0 == GrB_free(&level_do));
    TEST_CHECK(0 == GrB_free(&level_v));

    teardown();
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
    {NULL, NULL}
};
