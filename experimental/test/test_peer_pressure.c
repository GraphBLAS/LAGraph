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

const double coverage [] = {
    1.000000,
    0.653359,
    0.863014,
    0.544871,
    0.243590,
    0.833333
} ;

const double performance [] = {
    0.714286,
    0.989642,
    0.284034,
    0.974517,
    0.887701,
    0.866667
} ;

const double modularity [] = {
    0.000000,
    0.641262,
    0.084537,
    0.495322,
    0.158120,
    0.500000
} ;

//****************************************************************************
void test_peer_pressure(void)
{
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
        GxB_print (A, 5) ;

        // construct a directed graph G with adjacency matrix A
        OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        TEST_CHECK(A == NULL);

        // compute AT
        OK(LAGraph_Cached_AT(G, msg));

        GrB_Vector c = NULL;

        // compute clustering
        double cov, perf, mod ;
        OK(LAGr_PeerPressureClustering(&c, true, false, 0.0001, 50, G, msg));
        OK(LAGr_PartitionQuality(&cov, &perf, c, G->A, msg));
        OK(LAGr_Modularity(&mod, (double)1, c, G->A, msg));

        // GrB_Index n;
        // OK(GrB_Vector_size(&n, c));
        // LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT;

        // printf("\ncdlp (known result):\n");
        // OK(LAGraph_Vector_Print(cgood, pr, stdout, msg));
        // bool ok = false;
        // OK(LAGraph_Vector_IsEqual(&ok, c, cgood, msg));
        // TEST_CHECK(ok);

        // printf("\peer pressure:\n");
        // OK(LAGraph_Vector_Print(c, pr, stdout, msg));
        bool ok_cov = false, ok_perf = false, ok_mod = false;
        printf ("coverage:   %g %g\n", cov, coverage[k]) ;
        printf ("perf:       %g %g\n", perf, performance[k]) ;
        printf ("modularity: %g %g\n", mod, modularity[k]) ;
        ok_cov = (fabs(cov - coverage[k]) < 1e-4) ? true : ok_cov;
        ok_perf = (fabs(perf - performance[k]) < 1e-4) ? true : ok_perf;
        ok_mod = (fabs(mod - modularity[k]) < 1e-4) ? true : ok_mod;

        TEST_CHECK(ok_cov);
        TEST_CHECK(ok_perf);
        TEST_CHECK(ok_mod);

        OK(GrB_free(&c));

        OK(LAGraph_Delete(&G, msg));
    }

    LAGraph_Finalize(msg);
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors(void)
{
    LAGraph_Init(msg);

    snprintf(filename, LEN, LG_DATA_DIR "%s", "karate.mtx");
    FILE *f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    TEST_MSG("Loading of adjacency matrix failed");

    // construct an undirected graph G with adjacency matrix A
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL);

    GrB_Vector c = NULL;
    bool normalize = false, make_undirected = false;
    double thresh = 1e-5;
    int max_iter = 100;

    // c is NULL
    GrB_Info result = LAGr_PeerPressureClustering(
        &c, normalize, make_undirected, thresh, max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_NULL_POINTER);

    // G has no AT
    G->AT = NULL;
    make_undirected = true;
    G->kind = LAGraph_ADJACENCY_DIRECTED;
    G->is_symmetric_structure = LAGraph_FALSE;
    result = LAGr_PeerPressureClustering(&c, normalize, make_undirected,
                                            thresh, max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_NOT_CACHED);

    OK(LAGraph_Delete(&G, msg));
    LAGraph_Finalize(msg);
}

//****************************************************************************

TEST_LIST = {{"peer_pressure", test_peer_pressure},
             {"peer_pressure_errors", test_errors},
             {NULL, NULL}};
