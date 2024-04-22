//------------------------------------------------------------------------------
// LAGr_PartitionQuality: computes coverage and performance of a clustering
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

// The coverage of a graph clustering C ( Cov(C) ) is defined as the ratio of
// intra-cluster edges to the total edges in a graph. The performance of a graph
// clustering C ( Perf(C) ) is defined as the ratio of intra-cluster edges and
// inter-cluster non-edges to the total number of edges in a graph. These are
// very simple cluster quality metrics which can be used to evaluate the quality
// of a clustering, primarily based on the idea of intra-cluster density.

// https://arxiv.org/abs/0906.0612 pp. 15

#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&trace);                                                      \
        GrB_free(&k);                                                          \
        GrB_free(&C);                                                          \
        GrB_free(&CA);                                                         \
        GrB_free(&AT);                                                         \
        GrB_free(&ONE_BOOL);                                                   \
    }

#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_PartitionQuality(
    // Outputs
    double *cov,  // coverage output, can be NULL
    double *perf, // performance output, can be NULL
    // Inputs
    GrB_Vector c, // input cluster vector
    GrB_Matrix A, // adjacency matrix of original graph from which the
                  // clustering was obtained
    char *msg)
{
    GrB_Vector trace = NULL;
    GrB_Vector k = NULL;
    GrB_Matrix C = NULL;
    GrB_Matrix CA = NULL;

    GrB_Scalar ONE_BOOL = NULL;

    GrB_Matrix AT = NULL;

    GrB_Index n, nedges;
    GrB_Index n_intraEdges, n_interEdges, n_interNonEdges, sum_k2;

    LG_ASSERT_MSG(cov != NULL || perf != NULL, GrB_NULL_POINTER,
                  "cov and perf "
                  "cannot both be NULL");

    // Delete self-edges, not relevant to these clustering metrics
    GRB_TRY(GrB_select(A, NULL, NULL, GrB_OFFDIAG, A, 0, NULL));

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_nvals(&nedges, A));

    // USE int64:
    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&CA, GrB_INT64, n, n));
    GRB_TRY(GrB_Vector_new(&trace, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&k, GrB_INT64, n));
    GRB_TRY(GrB_Scalar_new(&ONE_BOOL, GrB_BOOL));

    GRB_TRY(GrB_Scalar_setElement_BOOL(ONE_BOOL, (bool)1));

    // convert the cluster vector to a boolean matrix C where
    // C(i, j) = 1 if and only if vertex j is in cluster i
    GrB_Index *cI, *cX;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&cI, n, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&cX, n, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(cI, cX, &n, c));
    GRB_TRY(GxB_Matrix_build_Scalar(C, cX, cI, ONE_BOOL, n));
    GrB_Matrix_wait(C, GrB_MATERIALIZE);
    LAGraph_Free((void **)&cI, NULL);
    LAGraph_Free((void **)&cX, NULL);

    // FIXME: pass in G and do
    // bool_is_undirected = (G->kind == LAGraph_ADJACENCY_UNDIRECTED) ;

    // DON'T:
    // check if graph is undirected (that is, is A symmetric)
    bool is_undirected;
    GRB_TRY(GrB_Matrix_new(&AT, GrB_BOOL, n, n));
    GRB_TRY(GrB_transpose(AT, NULL, NULL, A, NULL));
    LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&is_undirected, A, AT, msg));

    // k = sum(C) .^ 2
    GRB_TRY(GrB_reduce(k, NULL, NULL, GrB_PLUS_MONOID_INT64, C, NULL));
    GRB_TRY(GrB_Vector_apply_BinaryOp2nd_INT64(k, NULL, NULL, GxB_POW_INT64, k,
                                               2, NULL));
    // sum_k2 = total number of possible intra-cluster edges
    GRB_TRY(GrB_reduce(&sum_k2, NULL, GrB_PLUS_MONOID_INT64, k, NULL));

    // Calculate actual number of intra-cluster edges
    GRB_TRY(GrB_mxm(CA, NULL, NULL, LAGraph_plus_one_int64, C, A, NULL));
    GRB_TRY(GrB_mxm(CA, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, CA, C,
                    GrB_DESC_RT1));
    GRB_TRY(GxB_Vector_diag(trace, CA, 0, NULL));

    GRB_TRY(
        GrB_reduce(&n_intraEdges, NULL, GrB_PLUS_MONOID_INT64, trace, NULL));

    if (perf)
    {
        // undirected case
        if (is_undirected)
        {
            nedges /= 2;
            n_intraEdges /= 2;
            n_interEdges = nedges - n_intraEdges;
            // All possible edges - intra cluster non-edges gives space of
            // possible inter-cluster edges. Then subtract the actual
            // inter-cluster edges to get the number of inter-cluster non-edges.
            n_interNonEdges =
                (n * (n - 1) / 2) - ((sum_k2 - n) / 2) - n_interEdges;
            (*perf) =
                (n_intraEdges + n_interNonEdges) * 1.0 / (n * (n - 1) / 2);
        }
        // directed case
        else
        {
            n_interEdges = nedges - n_intraEdges;
            n_interNonEdges = n * (n - 1) - (sum_k2 - n) - n_interEdges;
            (*perf) = (n_intraEdges + n_interNonEdges) * 1.0 / (n * (n - 1));
        }
    }

    if (cov)
    {
        (*cov) = n_intraEdges * 1.0 / nedges;
    }

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}
