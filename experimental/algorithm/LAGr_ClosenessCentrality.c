//------------------------------------------------------------------------------
// LAGr_ClosenessCentrality: Closeness centrality algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Karan Bhalla and Timothy A. Davis, Texas A&M University;
// Adapted and revised from GraphBLAS C API Spec, Appendix B.4.

//------------------------------------------------------------------------------

#define LG_FREE_WORK                                                   \
    {                                                                  \
        GrB_free(&distance_vector);                                    \
        if (incoming_adjacency != NULL &&                              \
            incoming_adjacency != G->A && incoming_adjacency != G->AT) \
        {                                                              \
            GrB_free(&incoming_adjacency);                             \
        }                                                              \
        LAGraph_Free((void **)&source_indices, NULL);                  \
    }

#define LG_FREE_ALL                   \
    {                                 \
        LG_FREE_WORK;                 \
        GrB_free(&centrality_vector); \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_ClosenessCentrality(
    // output:
    GrB_Vector *centrality,
    // input:
    LAGraph_Graph G,
    GrB_Vector sources, // target vertices to score; NULL/empty => all
    const GrB_Matrix D, // optional APSP matrix, D(i,j)=dist(i->j), or NULL
    char *msg)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    GrB_Info info;

    // centrality_vector(v) will store the final closeness score for v.
    GrB_Vector centrality_vector = NULL;

    // distance_vector(u) stores shortest path length from the current target
    // node to u in the incoming-edge view of the graph.
    GrB_Vector distance_vector = NULL;

    // incoming_adjacency is A' for directed graphs and A for undirected ones.
    GrB_Matrix incoming_adjacency = NULL;

    GrB_Index n = 0;
    GrB_Index source_count = 0;

    // If a subset of nodes is requested, source_indices holds those node ids.
    // We treat 'sources' as a set/mask vector and use tuple indices as IDs.
    GrB_Index *source_indices = NULL;

    LG_ASSERT(centrality != NULL, GrB_NULL_POINTER);
    (*centrality) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GRB_TRY(GrB_Matrix_nrows(&n, G->A));

    if (D != NULL)
    {
        GrB_Index nrows, ncols;
        GRB_TRY(GrB_Matrix_nrows(&nrows, D));
        GRB_TRY(GrB_Matrix_ncols(&ncols, D));
        LG_ASSERT(nrows == n && ncols == n, GrB_DIMENSION_MISMATCH);
    }

    bool use_all_nodes = (sources == NULL);
    if (!use_all_nodes)
    {
        GRB_TRY(GrB_Vector_nvals(&source_count, sources));
        use_all_nodes = (source_count == 0);
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Vector_new(&centrality_vector, GrB_FP64, n));
    if (use_all_nodes)
    {
        // dense output: every node gets an entry
        GRB_TRY(GrB_assign(centrality_vector, NULL, NULL,
                           0.0, GrB_ALL, n, NULL));
        source_count = n;
    }

    // Uuses Bellman-Ford 
    // Does not consume an externally-provided all-pairs shortest path matrix yet.
    LG_ASSERT_MSG(D == NULL, GrB_NOT_IMPLEMENTED,
                  "D input is not supported in this simplified version");

    // incoming-edge closeness: for directed graphs we need distances TO each
    // node, so we compute shortest paths on A'.
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        incoming_adjacency = G->A;
    }
    else if (G->AT != NULL)
    {
        incoming_adjacency = G->AT;
    }
    else
    {
        GrB_Type edge_type;
        GRB_TRY(GxB_Matrix_type(&edge_type, G->A));
        GRB_TRY(GrB_Matrix_new(&incoming_adjacency, edge_type, n, n));
        GRB_TRY(GrB_transpose(incoming_adjacency, NULL, NULL, G->A, NULL));
    }

    // Build list of requested target nodes from tuple indices in 'sources'.
    if (!use_all_nodes)
    {
        LG_TRY(LAGraph_Malloc((void **)&source_indices,
                              source_count, sizeof(GrB_Index), msg));
        GRB_TRY(GrB_Vector_extractTuples(source_indices, NULL,
                                         &source_count, sources));
        for (GrB_Index k = 0; k < source_count; k++)
        {
            LG_ASSERT(source_indices[k] < n, GrB_INVALID_INDEX);
        }
    }

    //--------------------------------------------------------------------------
    // compute centrality values
    //--------------------------------------------------------------------------

    for (GrB_Index k = 0; k < source_count; k++)
    {
        // node_to_score is the vertex whose incoming closeness we compute.
        GrB_Index node_to_score = use_all_nodes ? k : source_indices[k];

        // Compute shortest distances from node_to_score in the incoming-edge
        // graph.  This equals distances TO node_to_score in the original graph.
        info = LAGraph_BF_basic(&distance_vector, incoming_adjacency,
                                node_to_score);
        if (info == GrB_NO_VALUE)
        {
            LG_FREE_ALL;
            return (info);
        }
        if (info < GrB_SUCCESS)
        {
            GRB_CATCH(info);
        }

        // Do not include distance(node_to_score, node_to_score)=0 in the sum.
        GRB_TRY(GrB_Vector_removeElement(distance_vector, node_to_score));

        // Count how many nodes can reach node_to_score.
        GrB_Index reachable_count = 0;
        GRB_TRY(GrB_Vector_nvals(&reachable_count, distance_vector));

        if (reachable_count > 0)
        {
            // Sum all finite shortest-path distances to node_to_score.
            double distance_sum = 0;
            GRB_TRY(GrB_reduce(&distance_sum, NULL, GrB_PLUS_MONOID_FP64,
                               distance_vector, NULL));
            if (distance_sum > 0)
            {
                // Closeness = (# reachable nodes) / (sum of distances).
                double closeness_value = ((double)reachable_count) / distance_sum;
                GRB_TRY(GrB_Vector_setElement(centrality_vector,
                                              closeness_value, node_to_score));
            }
        }

        GRB_TRY(GrB_free(&distance_vector));
    }

    (*centrality) = centrality_vector;
    LG_FREE_WORK;

    return GrB_SUCCESS;
}
