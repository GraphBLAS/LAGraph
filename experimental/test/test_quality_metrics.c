//----------------------------------------------------------------------------
// LAGraph/src/test/test_quality_metrics.c: test cases for both Partition
// Quality and Modularity graph clustering metrics
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
GrB_Vector c = NULL;
#define LEN 512
char filename[LEN + 1];
char cluster_filename[LEN + 1];

typedef struct
{
    const char *name;
    const char *cluster_name;
} matrix_info;

const matrix_info files[] = {
    {"A.mtx", "A_cluster.mtx"},
    {"jagmesh7.mtx", "jagmesh7_cluster.mtx"},
    {"west0067.mtx", "west0067_cluster.mtx"}, // unsymmetric
    {"bcsstk13.mtx", "bcsstk13_cluster.mtx"},
    {"karate.mtx", "karate_cluster.mtx"},
    {"mcl.mtx", "mcl_cluster.mtx"},
    {"", ""},
};

const int nfiles = 6;

const double coverage[] = {1.000000, 0.653359, 0.181507,
                           0.048510, 0.243590, 0.833333};

const double performance[] = {0.714286, 0.989642, 0.841701,
                              0.977048, 0.887701, 0.866667};

const double modularity[] = {0.000000, 0.641262, 0.043324,
                             0.042696, 0.158120, 0.500000};

//****************************************************************************
void test_quality_metrics(void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Init(msg);

    for (int k = 0;; k++)
    {
        // load the matrix as A
        const char *aname = files[k].name;
        const char *aname_cluster = files[k].cluster_name;
        if (strlen(aname) == 0 || strlen(aname_cluster) == 0)
            break;
        printf("\n================================== %s:, %s\n", aname,
               aname_cluster);
        TEST_CASE(aname);
        TEST_CASE(aname_cluster);
        snprintf(filename, LEN, LG_DATA_DIR "%s", aname);
        FILE *f1 = fopen(filename, "r");
        TEST_CHECK(f1 != NULL);
        OK(LAGraph_MMRead(&A, f1, msg));
        fclose (f1) ;

        snprintf(cluster_filename, LEN, LG_DATA_DIR "%s", aname_cluster);
        FILE *f2 = fopen(cluster_filename, "r");
        TEST_CHECK(f2 != NULL);
        OK(LAGraph_MMRead((GrB_Matrix *)&c, f2, msg));
        fclose (f2) ;

        // construct a directed graph G with adjacency matrix A
        OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        TEST_CHECK(A == NULL);

        // compute is_symmetric_structure
        OK(LAGraph_Cached_IsSymmetricStructure(G, msg));
        TEST_CHECK(G->is_symmetric_structure != LAGRAPH_UNKNOWN);

        // compute quality metrics (coverage and performance)
        double cov, perf, mod;
        OK(LAGr_PartitionQuality(&cov, NULL, c, G, msg));
        OK(LAGr_PartitionQuality(NULL, &perf, c, G, msg));
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

        OK(LAGraph_Delete(&G, msg));
    }

    LAGraph_Finalize(msg);
#endif
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_partition_quality_errors(void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Init(msg);

    snprintf(filename, LEN, LG_DATA_DIR "%s", "west0067.mtx");
    FILE *f1 = fopen(filename, "r");
    TEST_CHECK(f1 != NULL);
    OK(LAGraph_MMRead(&A, f1, msg));
    TEST_MSG("Loading of adjacency matrix failed");
    fclose (f1) ;

    snprintf(cluster_filename, LEN, LG_DATA_DIR "%s", "west0067_cluster.mtx");
    FILE *f2 = fopen(cluster_filename, "r");
    TEST_CHECK(f2 != NULL);
    OK(LAGraph_MMRead((GrB_Matrix *)&c, f2, msg));
    TEST_MSG("Loading of cluster vector failed");
    fclose (f2) ;

    // construct an undirected graph G with adjacency matrix A
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL);

    double cov, perf;

    // both cov and perf are null
    GrB_Info result = LAGr_PartitionQuality(NULL, NULL, c, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_NULL_POINTER);

    // G->is_symmetric_structure is not cached
    G->is_symmetric_structure = LAGRAPH_UNKNOWN;
    result = LAGr_PartitionQuality(&cov, &perf, c, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_NOT_CACHED);

    OK(LAGraph_Delete(&G, msg));
    TEST_CHECK(G == NULL);

    // bad graph, G->A is null
    OK(LAGraph_New(&G, NULL, LAGraph_ADJACENCY_UNDIRECTED, msg));
    result = LAGr_PartitionQuality(&cov, &perf, c, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_INVALID_GRAPH);

    GrB_free (&c) ;
    OK(LAGraph_Delete(&G, msg));
    TEST_CHECK(G == NULL);

    LAGraph_Finalize(msg);
#endif
}

void test_modularity_errors(void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Init(msg);

    snprintf(filename, LEN, LG_DATA_DIR "%s", "west0067.mtx");
    FILE *f1 = fopen(filename, "r");
    TEST_CHECK(f1 != NULL);
    OK(LAGraph_MMRead(&A, f1, msg));
    TEST_MSG("Loading of adjacency matrix failed");
    fclose (f1) ;

    snprintf(cluster_filename, LEN, LG_DATA_DIR "%s", "west0067_cluster.mtx");
    FILE *f2 = fopen(cluster_filename, "r");
    TEST_CHECK(f2 != NULL);
    OK(LAGraph_MMRead((GrB_Matrix *)&c, f2, msg));
    TEST_MSG("Loading of cluster vector failed");
    fclose (f2) ;

    // construct an undirected graph G with adjacency matrix A
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL);

    double mod;
    double resolution = 1;

    // mod is NULL
    GrB_Info result = LAGr_Modularity(NULL, resolution, c, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_NULL_POINTER);

    // resolution parameter is negative
    resolution = -1;
    result = LAGr_Modularity(&mod, resolution, c, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == GrB_INVALID_VALUE);
    resolution = 1;

    OK(LAGraph_Delete(&G, msg));
    TEST_CHECK(G == NULL);

    // bad graph, G->A is null
    OK(LAGraph_New(&G, NULL, LAGraph_ADJACENCY_UNDIRECTED, msg));
    result = LAGr_Modularity(&mod, resolution, c, G, msg);
    printf("\nresult: %d %s\n", result, msg);
    TEST_CHECK(result == LAGRAPH_INVALID_GRAPH);

    GrB_free (&c) ;
    OK(LAGraph_Delete(&G, msg));
    TEST_CHECK(G == NULL);

    LAGraph_Finalize(msg);
#endif
}

//****************************************************************************

TEST_LIST = {{"quality_metrics", test_quality_metrics},
             {"partition_quality_errors", test_partition_quality_errors},
             {"modularity_errors", test_modularity_errors},
             {NULL, NULL}};
