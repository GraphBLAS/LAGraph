//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/peer_pressure_demo.c:
// benchmark for peer pressure clustering algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2023 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Cameron Quilici, Texas A&M University

//------------------------------------------------------------------------------

// Usage:  peer_pressure_demo < matrixmarketfile.mtx
//         peer_pressure_demo matrixmarketfile.mtx
//         peer_pressure_demo matrixmarketfile.grb

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_Xtest.h"

#define LG_FREE_ALL                                                            \
    {                                                                          \
        LAGraph_Delete(&G, NULL);                                              \
        GrB_free(&c);                                                          \
        GrB_free(&C);                                                          \
        GrB_free(&vpc);                                                        \
        GrB_free(&vpc_sorted);                                                 \
        GrB_free(&TRUE_BOOL);                                                  \
        LAGraph_Free((void **)&cI, NULL);                                      \
        LAGraph_Free((void **)&cX, NULL);                                      \
    }

#define IO 0

int main(int argc, char **argv)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg[LAGRAPH_MSG_LEN];

    GrB_Vector c = NULL;
    GrB_Vector vpc = NULL, vpc_sorted = NULL;
    GrB_Matrix C = NULL;
    GrB_Scalar TRUE_BOOL = NULL;
    LAGraph_Graph G = NULL;

    GrB_Index *cI = NULL, *cX = NULL;

    // start GraphBLAS and LAGraph
    bool burble = false;
    demo_init(burble);

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv[1] : "stdin";
    LAGRAPH_TRY(
        readproblem(&G, NULL, false, false, false, NULL, false, argc, argv));

    GrB_Index n, nvals;
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));
    GRB_TRY(GrB_Matrix_nvals(&nvals, G->A));
    
    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Vector_new(&vpc, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&vpc_sorted, GrB_INT64, n));
    GRB_TRY(GrB_Scalar_new(&TRUE_BOOL, GrB_BOOL));

    GRB_TRY(GrB_Scalar_setElement_BOOL(TRUE_BOOL, (bool)1));

    //--------------------------------------------------------------------------
    // run peer pressure clustering algorithm
    //--------------------------------------------------------------------------

    // compute check result
    c = NULL;
    double tt = LAGraph_WallClockTime();
    LAGRAPH_TRY(
        LAGr_PeerPressureClustering(&c, true, false, 0.0001, 50, G, msg));
    tt = LAGraph_WallClockTime() - tt;
    printf("peer pressure run time %g sec\n", tt);

    // GxB_print(c, GxB_SHORT);

    double cov, perf, mod;
    tt = LAGraph_WallClockTime();
    LAGRAPH_TRY(LAGr_PartitionQuality(&cov, &perf, c, G, msg));
    tt = LAGraph_WallClockTime() - tt;
    printf("\npartition quality run time %g sec\n\tcoverage    = "
           "%f\n\tperformance = %f\n",
           tt, cov, perf);

    tt = LAGraph_WallClockTime();
    LAGRAPH_TRY(LAGr_Modularity(&mod, (double)1, c, G, msg));
    tt = LAGraph_WallClockTime() - tt;
    printf("modularity run time %g sec\n\tmodularity  = %f\n", tt, mod);

    GrB_Index nclusters;
    GRB_TRY(GrB_Vector_nvals(&nclusters, c));

    //--------------------------------------------------------------------------
    // calculate vertices per cluster by constructing cluster matrix
    //--------------------------------------------------------------------------

    LAGRAPH_TRY(
        LAGraph_Malloc((void **)&cI, nclusters, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(
        LAGraph_Malloc((void **)&cX, nclusters, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(cI, (int64_t *) cX, &nclusters, c));
    GRB_TRY(GxB_Matrix_build_Scalar(C, cX, cI, TRUE_BOOL, nclusters));

    GRB_TRY(GrB_reduce(vpc, NULL, NULL, GrB_PLUS_MONOID_INT64, C, NULL));
    GRB_TRY(GxB_Vector_sort(vpc_sorted, NULL, GrB_GT_FP64, vpc, NULL));

    GxB_print(vpc_sorted, GxB_SHORT);

#if IO
    FILE *f = fopen("./data/pp_out.mtx", "w");
    LAGRAPH_TRY(LAGraph_MMWrite((GrB_Matrix)c, f, NULL, msg));
    fclose(f);
#endif

    LG_FREE_ALL;
    LAGRAPH_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
