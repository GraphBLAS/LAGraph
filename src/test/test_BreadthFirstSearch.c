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

    //printf("parent BFS (!pushpull) source = 30\n");
    //LAGraph_Vector_print_type(parent, GrB_INT32, 3, stdout, msg);

    retval = LAGraph_BreadthFirstSearch(NULL, &parent_do, G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_parents30(parent_do));

    //printf("parent BFS (pushpull) source = 30\n");
    //LAGraph_Vector_print_type(parent_do, GrB_INT32, 3, stdout, msg);

    TEST_CHECK(0 == GrB_free(&parent));
    TEST_CHECK(0 == GrB_free(&parent_do));

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
    TEST_CHECK(0 == GrB_free(&level));

    retval = LAGraph_BreadthFirstSearch(&level_do, NULL, G, 30, true, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);
    TEST_CHECK(check_karate_levels30(level_do));
    TEST_CHECK(0 == GrB_free(&level_do));

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

    TEST_CHECK(0 == GrB_free(&parent_do));
    TEST_CHECK(0 == GrB_free(&level_do));

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
