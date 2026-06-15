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

// TODO: ready to consider for src
// TODO: document the input vector c that defines the cluster assignment

// The coverage of a graph clustering C ( Cov(C) ) is defined as the ratio of
// intra-cluster edges to the total edges in a graph. The performance of a graph
// clustering C ( Perf(C) ) is defined as the ratio of intra-cluster edges and
// inter-cluster non-edges to the total number of edges in a graph. These are
// very simple cluster quality metrics which can be used to evaluate the quality
// of a clustering, primarily based on the idea of intra-cluster density.

// Note that coverage and performance are counting problems. Therefore, if a
// weighted graph is passed to this function, edge weights are ignored.

// https://arxiv.org/abs/0906.0612 pp. 15

#include "GraphBLAS.h"
#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&k);                                                          \
        GrB_free(&x);                                                          \
        GrB_free(&C);                                                          \
        GrB_free(&ONE_INT64);                                                  \
    }

#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

// FIXME: make cov, perf outputs GrB_Scalar

int LAGr_PartitionQuality(
    // Outputs
    double *cov,  // coverage output, can be NULL
    double *perf, // performance output, can be NULL
    // Inputs
    GrB_Vector c,    // input cluster vector
    LAGraph_Graph G, // original graph from which the clustering was obtained
    char *msg)
{
#if LAGRAPH_SUITESPARSE
    GrB_Vector k = NULL;
    GrB_Vector x = NULL;
    GrB_Matrix C = NULL;

    GrB_Scalar ONE_INT64 = NULL;

    GrB_Index n, nedges;
    GrB_Index n_intraEdges, n_interEdges, n_interNonEdges, sum_k2;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_ASSERT_MSG(cov != NULL || perf != NULL, GrB_NULL_POINTER,
                  "cov and perf "
                  "cannot both be NULL");

    LG_TRY(LAGraph_CheckGraph(G, msg));

    LG_ASSERT_MSG (G->is_symmetric_structure != LAGRAPH_UNKNOWN,
            LAGRAPH_NOT_CACHED,
            "G->is_symmetric_structure is required") ;

    LG_ASSERT_MSG (G->nself_edges != LAGRAPH_UNKNOWN,
            LAGRAPH_NOT_CACHED,
            "G->nself_edges is required") ;

    GrB_Matrix A = G->A;


#if 0
    FILE *f = fopen("./data/pp_sanitized_data.mtx", "w");
    LAGRAPH_TRY(LAGraph_MMWrite(A, f, NULL, msg));
    fclose(f);
#endif

    GRB_TRY (GrB_Matrix_nrows (&n, A));
    GRB_TRY (GrB_Matrix_nvals (&nedges, A));

    GRB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n));
    GRB_TRY (GrB_Vector_new (&k, GrB_INT64, n));
    GRB_TRY (GrB_Vector_new (&x, GrB_BOOL, n));
    GRB_TRY (GrB_Scalar_new (&ONE_INT64, GrB_INT64));

    GRB_TRY (GrB_assign(x, NULL, NULL, (bool) true, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Scalar_setElement_BOOL(ONE_INT64, (int64_t)1));

    // convert the cluster vector to a boolean matrix C where
    // C(i, j) = 1 if and only if vertex j is in cluster i
    GrB_Index *cI, *cX;
    LAGRAPH_TRY (LAGraph_Malloc ((void **)&cI, n, sizeof(GrB_Index), msg));
    LAGRAPH_TRY (LAGraph_Malloc ((void **)&cX, n, sizeof(GrB_Index), msg));
    GRB_TRY (GrB_Vector_extractTuples_INT64 (cI, (int64_t *) cX, &n, c));
    GRB_TRY (GxB_Matrix_build_Scalar (C, cX, cI, ONE_INT64, n));
    LAGraph_Free((void **)&cI, NULL);
    LAGraph_Free((void **)&cX, NULL);

    bool is_undirected = (G->is_symmetric_structure == LAGraph_TRUE);

    // k = sum(C) .^ 2
    GRB_TRY (GrB_mxv (k, NULL, NULL, GxB_PLUS_PAIR_INT64, C, x, NULL));
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64(
        k, NULL, NULL, GxB_POW_INT64, k, 2, NULL));
    // sum_k2 = total number of possible intra-cluster edges
    GRB_TRY (GrB_reduce (&sum_k2, NULL, GrB_PLUS_MONOID_INT64, k, NULL));

    // Calculate actual number of intra-cluster edges, if A is weighted
    // then we ignore the weights and just count the number of edges as
    // performance and coverage are inherently counting problems

    // Count the number of edges per node in the cluster originating at a node
    // in the cluster
    GRB_TRY (GrB_mxm (C, C, NULL, GxB_PLUS_PAIR_INT64, C, A, GrB_DESC_S));

    // sum the intra-cluster edges per node for all nodes
    GRB_TRY (
        GrB_reduce (&n_intraEdges, NULL, GrB_PLUS_MONOID_INT64, C, NULL)) ;

    if (G->nself_edges > 0)
    {
        n_intraEdges -= G->nself_edges;
        nedges -= G->nself_edges;
    }

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
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
