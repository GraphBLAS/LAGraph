//------------------------------------------------------------------------------
// LAGr_Modularity: computes modularity of a graph clustering
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Cameron Quilici, Texas A&M University

//------------------------------------------------------------------------------

// TODO: ready to consider for src

// The modularity (Q) of a graph clustering C is defined as (directed case):
//
// Q = \sum{c = 1}^n [\frac{L_c}{m} - \gamma(\frac{k_c^{in} k_c^{out}}{m})^2]
//
// where the sum iterates over all communities c = 1,...,n, L_c is the the
// number of total edges in cluster c, m is the total number of edges in the
// graph, k_c^{in} and k_c^{out} are the total in and out degrees of cluster c,
// and \gamma is the resolution parameter which controls the relative importance
// of the intra-cluster and inter-cluster edges.

// Modularity works by comparing the intra-cluster density of a particular
// clustering to that of a random graph with the same degree distribution. The
// range of Q is [-.5, 1]. When Q ~ 0, this indicates that the clustering is
// not much better than at random whereas Q ~ 1 indicates a strong community
// structure.

// https://arxiv.org/abs/0906.0612 pp. 15-16

#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&l);                                                          \
        GrB_free(&vmask);                                                      \
        GrB_free(&k_in);                                                       \
        GrB_free(&k_out);                                                      \
        GrB_free(&in_degree);                                                  \
        GrB_free(&out_degree);                                                 \
        GrB_free(&A);                                                          \
        GrB_free(&CA);                                                         \
        GrB_free(&C);                                                          \
        GrB_free(&ONE_INT64);                                                  \
        LAGraph_Free((void **)&lX, NULL);                                      \
        LAGraph_Free((void **)&k_outX, NULL);                                  \
        LAGraph_Free((void **)&k_inX, NULL);                                   \
    }

#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_Modularity(
    // Outputs
    double *mod_handle, // output modularity
    // Inputs
    double resolution, // resolution parameter
    GrB_Vector c,      // input cluster vector
    LAGraph_Graph G,   // original graph from which clustering was obtained
    char *msg)
{
#if LAGRAPH_SUITESPARSE
    GrB_Vector l = NULL;
    GrB_Vector vmask = NULL;
    GrB_Vector k_in = NULL, k_out = NULL;
    GrB_Vector out_degree = NULL, in_degree = NULL;
    GrB_Matrix C = NULL, CA = NULL, A = NULL;
    GrB_Scalar ONE_INT64 = NULL;

    GrB_Index *lX = NULL, *k_outX = NULL, *k_inX = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_ASSERT_MSG(mod_handle != NULL, GrB_NULL_POINTER, "mod_handle is NULL");
    LG_ASSERT_MSG(resolution >= 0, GrB_INVALID_VALUE,
                  "resolution parameter must be non-negative");

    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Index n, nedges;
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));

    // Do not consider edge weights for modularity
    // FUTURE: There is a way to consider edge weights in modularity calculations,
    // so could give users the option to pass in a 'weights' parameter specifying 
    // that edge weights should be used in the calculations.
    GRB_TRY(GrB_Matrix_new(&A, GrB_INT64, n, n));
    GRB_TRY(GrB_apply(A, NULL, NULL, GxB_ONE_INT64, G->A, NULL));

    // remove self-edges, not relevant to clustering metrics
    GRB_TRY(GrB_select(A, NULL, NULL, GrB_OFFDIAG, A, 0, NULL));

    GRB_TRY(GrB_Matrix_nvals(&nedges, A));

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&C, GrB_INT64, n, n));
    GRB_TRY(GrB_Matrix_new(&CA, GrB_INT64, n, n));
    GRB_TRY(GrB_Vector_new(&l, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&vmask, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&k_in, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&k_out, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&out_degree, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&in_degree, GrB_INT64, n));
    GRB_TRY(GrB_Scalar_new(&ONE_INT64, GrB_INT64));

    GRB_TRY(GrB_Scalar_setElement_INT64(ONE_INT64, (int64_t)1));

    // Convert the cluster vector to a uint64_t matrix C where
    // C(i, j) = 1 if and only if vertex j is in cluster i
    GrB_Index *cI, *cX;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&cI, n, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&cX, n, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(cI, (int64_t *) cX, &n, c));
    GRB_TRY(GxB_Matrix_build_Scalar(C, cX, cI, ONE_INT64, n));
    GrB_Matrix_wait(C, GrB_MATERIALIZE);
    LAGraph_Free((void **)&cI, NULL);
    LAGraph_Free((void **)&cX, NULL);

    // Calculate actual number of intra-cluster edges
    GRB_TRY(GrB_mxm(CA, NULL, NULL, LAGraph_plus_one_int64, C, A, NULL));
    GRB_TRY(GrB_mxm(CA, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, CA, C,
                    GrB_DESC_T1));
    GRB_TRY(GxB_Vector_diag(l, CA, 0, NULL));

    // Calculate the combined degree for each cluster
    GRB_TRY(GrB_reduce(out_degree, NULL, NULL, GrB_PLUS_MONOID_INT64, A, NULL));
    GRB_TRY(GrB_reduce(in_degree, NULL, NULL, GrB_PLUS_MONOID_INT64, A,
                       GrB_DESC_T0));
    GRB_TRY(GrB_mxv(k_out, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C,
                    out_degree, NULL));
    GRB_TRY(GrB_mxv(k_in, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C,
                    in_degree, NULL));

    // vmask (i) = 0 if cluster i is non-empty (has any vertices)
    GRB_TRY(GrB_reduce(vmask, NULL, NULL, GrB_LOR_MONOID_BOOL, C, NULL));
    GRB_TRY(GrB_apply(vmask, vmask, NULL, GxB_LNOT_BOOL, vmask, NULL));

    // If any of the above vectors have fewer entries than nclusters, this means
    // that there are singleton clusters with one vertex/no out-degree/no
    // in-degree. So we need to add explicit zeros wherever values are missing
    // for further calculations.
    GrB_Index nclusters, nl, nk_out, nk_in;
    GRB_TRY(GrB_Vector_nvals(&nclusters, vmask));
    GRB_TRY(GrB_Vector_nvals(&nl, l));
    GRB_TRY(GrB_Vector_nvals(&nk_out, l));
    GRB_TRY(GrB_Vector_nvals(&nk_in, l));

    if (nclusters != nl)
    {
        GRB_TRY(GrB_assign(l, l, NULL, vmask, GrB_ALL, nclusters, GrB_DESC_SC));
    }
    if (nclusters != nk_out)
    {
        GRB_TRY(GrB_assign(k_out, k_out, NULL, vmask, GrB_ALL, nclusters,
                           GrB_DESC_SC));
    }
    if (nclusters != nk_in)
    {
        GRB_TRY(GrB_assign(k_in, k_in, NULL, vmask, GrB_ALL, nclusters,
                           GrB_DESC_SC));
    }

    // Extract actual values of l, k_out, and k_in for modularity calculations
    LAGRAPH_TRY(
        LAGraph_Malloc((void **)&lX, nclusters, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(
        LAGraph_Malloc((void **)&k_outX, nclusters, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(
        LAGraph_Malloc((void **)&k_inX, nclusters, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, (int64_t *) lX, &nclusters, l));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, (int64_t *) k_outX, &nclusters, k_out));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, (int64_t *) k_inX, &nclusters, k_in));

    GrB_Index m, out_degree_sum, in_degree_sum, L_c;
    GRB_TRY(GrB_reduce(&out_degree_sum, NULL, GrB_PLUS_MONOID_INT64, out_degree,
                       NULL));

    m = out_degree_sum;
    double norm = 1.0 / (m * m);

    // compute modularity
    double mod = 0.0;
    for (int c = 0; c < nclusters; c++)
    {
        mod += (1.0 * lX[c] / nedges) -
               (resolution * ((k_outX[c] * k_inX[c]) * norm));
    }

    (*mod_handle) = mod;
    LG_FREE_WORK;
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
