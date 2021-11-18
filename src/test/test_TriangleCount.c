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

#include <LAGraph_test.h>
#include <graph_zachary_karate.h>

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    uint64_t ntriangles ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    {     45, "karate.mtx" }, 
    {     11, "A.mtx" }, 
    {   2016, "jagmesh7.mtx" }, 
    {      6, "ldbc-cdlp-undirected-example.mtx" }, 
    {      4, "ldbc-undirected-example.mtx" }, 
    {      5, "ldbc-wcc-example.mtx" }, 
    {      0, "LFAT5.mtx" }, 
    { 342300, "bcsstk13.mtx" }, 
    {      0, "tree-example.mtx" }, 
    {      0, "" }, 
} ;

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

    retval = LAGraph_Property_NDiag(G, msg);
    TEST_CHECK(retval == 0);
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK(G->ndiag == 0);
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
    TEST_CHECK(retval == -106);  // should fail (rowdegrees needs to be defined)
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
    TEST_CHECK(retval == -106);  // should fail (rowdegrees needs to be defined)
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
    TEST_CHECK(retval == -106);  // should fail (rowdegrees needs to be defined)
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
    TEST_CHECK(retval == -106);  // should fail (rowdegrees needs to be defined)
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
void test_TriangleCount(void)
{
    setup();

    uint64_t ntriangles = 0UL;
    int retval = LAGraph_TriangleCount(&ntriangles, G, msg);
    TEST_CHECK(retval == 0);  // should not fail (rowdegrees will be calculated)
    TEST_MSG("retval = %d (%s)", retval, msg);

    TEST_CHECK( ntriangles == 45 );
    TEST_MSG("numtri = %ld", ntriangles);

    #if LG_SUITESPARSE
    OK (LG_check_tri (&ntriangles, G, msg)) ;
    TEST_CHECK( ntriangles == 45 );
    #endif

    teardown();
}

//****************************************************************************
void test_TriangleCount_brutal (void)
{
    LAGraph_Init(msg);
    GrB_Matrix A = NULL ;
    GrB_Type atype = NULL ;
    #if LG_SUITESPARSE
    printf ("\nResults compared with LG_check_tri:\n") ;
    #else
    printf ("\n") ;
    #endif

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        uint64_t ntriangles = files [k].ntriangles ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // delete any diagonal entries
        OK (LAGraph_DeleteDiag (G, msg)) ;
        TEST_CHECK (G->ndiag == 0) ;
        OK (LAGraph_DeleteDiag (G, msg)) ;
        TEST_CHECK (G->ndiag == 0) ;

        // get the # of triangles
        uint64_t nt0, nt1 ;
        OK (LAGraph_TriangleCount (&nt1, G, msg)) ;
        printf ("# triangles: %6lu Matrix: %s\n", nt1, aname) ;
        TEST_CHECK (nt1 == ntriangles) ;
        #if LG_SUITESPARSE
        OK (LG_check_tri (&nt0, G, msg)) ;
        TEST_CHECK (nt0 == nt1) ;
        #endif

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize(msg);
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
    {"TriangleCount", test_TriangleCount},
    {"TriangleCount_brutal", test_TriangleCount_brutal},
    {NULL, NULL}
};
