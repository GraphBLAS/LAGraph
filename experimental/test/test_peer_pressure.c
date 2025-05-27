//----------------------------------------------------------------------------
// LAGraph/src/test/test_peer_pressure.c: test cases for Peer Pressure
// ----------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Cameron Quilici, Texas A&M University

//-----------------------------------------------------------------------------

#include <acutest.h>
#include <stdio.h>

#include "LG_Xtest.h"
#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg[LAGRAPH_MSG_LEN];

LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
#define LEN 512
char filename[LEN + 1];

typedef struct
{
    const char *name;
} matrix_info;

const matrix_info files[] = {
    {"A.mtx"},        {"jagmesh7.mtx"}, {"west0067.mtx"}, // unsymmetric
    {"bcsstk13.mtx"}, {"karate.mtx"},   {"mcl.mtx"},      {""},
};

const int nfiles = 6;

const double coverage[] = {
    1.000000, 0.653359, 0.181507, 0.048510, 0.243590, 0.833333,
    // Start config 2
    1.000000, 0.644804, 0.123288, 0.695750, 1.000000, 0.722222};

const double performance[] = {
    0.714286, 0.989642, 0.841701, 0.977048, 0.887701, 0.866667,
    // Start config 2
    0.714286, 0.992349, 0.914518, 0.934843, 0.139037, 0.777778};

const double modularity[] = {
    0.000000, 0.641262, 0.043324, 0.042696, 0.158120, 0.500000,
    // Start config 2
    0.000000, 0.634677, 0.078228, 0.596324, 0.000000, 0.351852};

//****************************************************************************
void test_peer_pressure(void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Init(msg);

    for (int k = 0;; k++)
    {
        // load the matrix as A
        const char *aname = files[k].name;
        if (strlen(aname) == 0)
            break;
        printf("\n================================== %s:\n", aname);
        TEST_CASE(aname);
        snprintf(filename, LEN, LG_DATA_DIR "%s", aname);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        // GxB_print (A, 5) ;
        fclose (f) ;

        // construct a directed graph G with adjacency matrix A
        OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        TEST_CHECK(A == NULL);

        // compute AT
        OK(LAGraph_Cached_AT(G, msg));
        // Needed to compute quality metrics
        OK(LAGraph_Cached_IsSymmetricStructure(G, msg));

        GrB_Vector c = NULL;

        // compute clustering
        double cov, perf, mod;
        OK(LAGr_PeerPressureClustering(&c, true, false, 0.0001, 50, G, msg));
        OK(LAGr_PartitionQuality(&cov, &perf, c, G, msg));
        OK(LAGr_Modularity(&mod, (double)1, c, G, msg));
        GrB_free (&c) ;

        bool ok_cov = false, ok_perf = false, ok_mod = false;
        printf("\nConfiguration 1:\n");
        printf("coverage:   %g %g\n", cov, coverage[k]);
        printf("perf:       %g %g\n", perf, performance[k]);
        printf("modularity: %g %g\n", mod, modularity[k]);
        ok_cov = (fabs(cov - coverage[k]) < 1e-4) ? true : ok_cov;
        ok_perf = (fabs(perf - performance[k]) < 1e-4) ? true : ok_perf;
        ok_mod = (fabs(mod - modularity[k]) < 1e-4) ? true : ok_mod;

        TEST_CHECK(ok_cov);
        TEST_CHECK(ok_perf);
        TEST_CHECK(ok_mod);

        OK(LAGr_PeerPressureClustering(&c, false, true, 0.0001, 50, G, msg));
        OK(LAGr_PartitionQuality(&cov, &perf, c, G, msg));
        OK(LAGr_Modularity(&mod, (double)1, c, G, msg));

        ok_cov = false, ok_perf = false, ok_mod = false;
        printf("\nConfiguration 2:\n");
        printf("coverage:   %g %g\n", cov, coverage[k + nfiles]);
        printf("perf:       %g %g\n", perf, performance[k + nfiles]);
        printf("modularity: %g %g\n", mod, modularity[k + nfiles]);
        ok_cov = (fabs(cov - coverage[k + nfiles]) < 1e-4) ? true : ok_cov;
        ok_perf =
            (fabs(perf - performance[k + nfiles]) < 1e-4) ? true : ok_perf;
        ok_mod = (fabs(mod - modularity[k + nfiles]) < 1e-4) ? true : ok_mod;

        TEST_CHECK(ok_cov);
        TEST_CHECK(ok_perf);
        TEST_CHECK(ok_mod);

        OK(GrB_free(&c));

        OK(LAGraph_Delete(&G, msg));
    }

    LAGraph_Finalize(msg);
#endif
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors(void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Init(msg);

    snprintf(filename, LEN, LG_DATA_DIR "%s", "karate.mtx");
    FILE *f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    TEST_MSG("Loading of adjacency matrix failed");
    fclose (f) ;

    // construct an undirected graph G with adjacency matrix A
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL);

    GrB_Vector c = NULL;
    bool normalize = false, make_undirected = false;
    double thresh = 1e-5;
    int max_iter = 100;

    // c is NULL
    GrB_Info result = LAGr_PeerPressureClustering(
        NULL, normalize, make_undirected, thresh, max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_NULL_POINTER);

    // G has no AT
    G->AT = NULL;
    make_undirected = true;
    G->kind = LAGraph_ADJACENCY_DIRECTED;
    G->is_symmetric_structure = LAGraph_FALSE;
    result = LAGr_PeerPressureClustering(&c, normalize, make_undirected, thresh,
                                         max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_NOT_CACHED);
    GrB_free (&c) ;

    G->kind = LAGraph_ADJACENCY_UNDIRECTED;
    G->is_symmetric_structure = LAGraph_FALSE;
    result = LAGr_PeerPressureClustering(&c, normalize, make_undirected, thresh,
                                         max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_NOT_CACHED);
    GrB_free (&c) ;

    OK(LAGraph_Delete(&G, msg));
    LAGraph_Finalize(msg);
#endif
}

//****************************************************************************

TEST_LIST = {{"peer_pressure", test_peer_pressure},
             {"peer_pressure_errors", test_errors},
             {NULL, NULL}};
