//----------------------------------------------------------------------------
// LAGraph/src/test/test_TriangleCount.cpp: test cases for triangle
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

#include <LAGraph.h>
#include <graph_zachary_karate.h>

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

//****************************************************************************
void setup(void)
{
    LAGraph_Init(msg);
    int retval;

    GrB_Matrix A = NULL;

    GrB_Matrix_new(&A, GrB_UINT32, ZACHARY_NUM_NODES, ZACHARY_NUM_NODES);
    GrB_Matrix_build(A, ZACHARY_I, ZACHARY_J, ZACHARY_V, ZACHARY_NUM_EDGES,
                     GrB_LOR);

    retval = LAGraph_New(&G, &A, GrB_UINT32, LAGRAPH_ADJACENCY_UNDIRECTED, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    G->ndiag = 0;         // just trust that the graph has no self loops.
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
void test_TriangleCount_Methods1(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    // no presort
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 1, NULL, msg);
    fprintf(stderr, "ret, n = %d %ld \n", retval, ntriangles);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    // with presort
    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 1, &presort, msg);

    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_Methods2(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 2, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 2, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_Methods3(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 3, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_Methods4(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 4, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );

#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 4, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_Methods5(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 5, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 5, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_Methods6(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 6, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 6, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_Methods(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
//****************************************************************************


//****************************************************************************
void test_TriangleCount_vanilla1(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    // no presort
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 1, NULL, msg);
    fprintf(stderr, "ret, n = %d %ld \n", retval, ntriangles);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    // with presort
    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 1, &presort, msg);

    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_vanilla2(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 2, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 2, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_vanilla3(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 3, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_vanilla4(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 4, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 4, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_vanilla5(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 5, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 5, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
void test_TriangleCount_vanilla6(void)
{
    setup();
    int retval;
    uint64_t ntriangles = 0UL;

#if 0
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 6, NULL, msg);
    TEST_CHECK(retval == 0);
    TEST_CHECK( ntriangles == 45 );
#endif

    int presort = 2;
    ntriangles = 0UL;
    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 6, &presort, msg);
    TEST_CHECK(retval == -6);  // should fail (rowdegrees needs to be defined)
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_Property_RowDegree(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    retval = LAGraph_TriangleCount_vanilla(&ntriangles, G, 3, &presort, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    teardown();
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"TriangleCount_Methods1", test_TriangleCount_Methods1},
    {"TriangleCount_Methods2", test_TriangleCount_Methods2},
    {"TriangleCount_Methods3", test_TriangleCount_Methods3},
    {"TriangleCount_Methods4", test_TriangleCount_Methods4},
    {"TriangleCount_Methods5", test_TriangleCount_Methods5},
    {"TriangleCount_Methods6", test_TriangleCount_Methods6},
    {"TriangleCount_vanilla1", test_TriangleCount_vanilla1},
    {"TriangleCount_vanilla2", test_TriangleCount_vanilla2},
    {"TriangleCount_vanilla3", test_TriangleCount_vanilla3},
    {"TriangleCount_vanilla4", test_TriangleCount_vanilla4},
    {"TriangleCount_vanilla5", test_TriangleCount_vanilla5},
    {"TriangleCount_vanilla6", test_TriangleCount_vanilla6},
    {NULL, NULL}
};
