//----------------------------------------------------------------------------
// LAGraph/src/test/test_mcl.c: test cases for Markov Clustering
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

const double coverage[] = {1.000000, 0.635932, 0.784247,
                           0.089424, 0.871795, 0.888889};

const double performance[] = {0.714286, 0.990614, 0.282678,
                              0.975945, 0.611408, 0.622222};

const double modularity[] = {0.000000, 0.624182, 0.033355,
                             0.083733, 0.359961, 0.339506};

//****************************************************************************
void test_mcl(void)
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
        OK(LAGr_MarkovClustering(&c, 2, 2, 0.0001, 1e-8, 100, G, msg));
        OK(LAGr_PartitionQuality(&cov, &perf, c, G, msg));
        OK(LAGr_Modularity(&mod, (double)1, c, G, msg));

        bool ok_cov = false, ok_perf = false, ok_mod = false;
        printf("coverage:   %g %g\n", cov, coverage[k]);
        printf("perf:       %g %g\n", perf, performance[k]);
        printf("modularity: %g %g\n", mod, modularity[k]);
        ok_cov = (fabs(cov - coverage[k]) < 1e-4) ? true : ok_cov;
        ok_perf = (fabs(perf - performance[k]) < 1e-4) ? true : ok_perf;
        ok_mod = (fabs(mod - modularity[k]) < 1e-4) ? true : ok_mod;

        TEST_CHECK(ok_cov);
        TEST_CHECK(ok_perf);
        TEST_CHECK(ok_mod);
        OK(GrB_free(&c));

        if (k != 3)
        {
            // compute clustering with higher e parameter (expansion coef)
            printf ("\nWith e=4:\n") ;
            OK(LAGr_MarkovClustering(&c, 4, 2, 0.0001, 1e-8, 100, G, msg));
            OK(LAGr_PartitionQuality(&cov, &perf, c, G, msg));
            OK(LAGr_Modularity(&mod, (double)1, c, G, msg));
            printf("coverage:   %g\n", cov);
            printf("perf:       %g\n", perf);
            printf("modularity: %g\n", mod);
            OK(GrB_free(&c));

            // compute clustering with high pruning threshold
            printf ("\nWith high pruning threshold:\n") ;
            OK (LAGr_MarkovClustering(&c, 4, 2, 0.005, 1e-8, 100, G, msg));
            OK (LAGr_PartitionQuality(&cov, &perf, c, G, msg));
            OK (LAGr_Modularity(&mod, (double)1, c, G, msg));
            printf("coverage:   %g\n", cov);
            printf("perf:       %g\n", perf);
            printf("modularity: %g\n", mod);
            OK(GrB_free(&c));
        }

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
    int e = 2, i = 2, max_iter = 50;
    double prune_thresh = 0.0001, conv_thresh = 1e-8;

    // c is NULL
    GrB_Info result = LAGr_MarkovClustering(NULL, e, i, prune_thresh,
                                            conv_thresh, max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_NULL_POINTER);

    // e is less than 2
    e = -100;
    result = LAGr_MarkovClustering(&c, e, i, prune_thresh, conv_thresh,
                                   max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_INVALID_VALUE);

    OK(LAGraph_Delete(&G, msg));
    TEST_CHECK(G == NULL);

    // bad graph, G->A is null
    OK(LAGraph_New(&G, NULL, LAGraph_ADJACENCY_UNDIRECTED, msg));
    result = LAGr_MarkovClustering(&c, e, i, prune_thresh, conv_thresh,
                                   max_iter, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_INVALID_GRAPH);

    OK(LAGraph_Delete(&G, msg));
    TEST_CHECK(G == NULL);

    LAGraph_Finalize(msg);
#endif
}

//****************************************************************************

TEST_LIST = {{"mcl", test_mcl},
             {"mcl_errors", test_errors},
             {NULL, NULL}};
