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

// Contributed by Karan Bhalla and Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#define LG_FREE_WORK                                  \
    {                                                 \
        GrB_free(&distance_vector);                   \
        GrB_free(&bfs_level);                         \
        GrB_free(&inf_scalar);                        \
        GrB_free(&dist_sums);                         \
        GrB_free(&reachable_counts);                  \
        GrB_free(&D_APSP);                            \
        GrB_free(&A_work);                            \
        GrB_free(&x);                                 \
        GrB_free(&Delta);                             \
        LAGraph_Delete(&G_AT, NULL);                  \
        LAGraph_Free((void **)&source_indices, NULL); \
    }

#define LG_FREE_ALL                   \
    {                                 \
        LG_FREE_WORK;                 \
        GrB_free(&centrality_vector); \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

typedef enum
{
    CC_BFS,           // unweighted BFS (unit edge weights)
    CC_SSSP,          // delta-stepping SSSP (non-negative weights)
    CC_BELLMAN_FORD,  // Bellman-Ford (general weights, negative-cycle check)
    CC_FLOYD_WARSHALL // Floyd-Warshall (FW) APSP
} cc_algo_t;

//------------------------------------------------------------------------------
// LAGr_ClosenessCentrality
//------------------------------------------------------------------------------

int LAGr_ClosenessCentrality(
    // output:
    GrB_Vector *centrality,
    // input:
    LAGraph_Graph G,
    GrB_Vector sources,      // nodes to score; NULL or empty => all nodes
    bool use_weights,        // if true, use edge weights in shortest paths
    bool use_floyd_warshall, // if true and sources==NULL, use APSP via FW
    char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;
    GrB_Info info;
    GrB_Vector centrality_vector = NULL;
    GrB_Vector distance_vector = NULL;  // per-node BF / SSSP result
    GrB_Vector bfs_level = NULL;        // BFS: level (distance) of each node
    GrB_Scalar inf_scalar = NULL;       // INFINITY threshold for SSSP pruning
    GrB_Vector dist_sums = NULL;        // FW: row sums of D
    GrB_Vector reachable_counts = NULL; // FW: row non-zero counts of D
    GrB_Matrix D_APSP = NULL;           // FW: APSP result
    GrB_Matrix A_work = NULL;           // AT copy for G_AT
    LAGraph_Graph G_AT = NULL;          // temporary graph wrapping AT (directed only)
    GrB_Vector x = NULL;                // FW: dense all-zeros vector for row counting
    GrB_Scalar Delta = NULL;            // delta-stepping step size
    GrB_Index *source_indices = NULL;
    GrB_Index n = 0;
    GrB_Index source_count = 0;

    LG_ASSERT(centrality != NULL, GrB_NULL_POINTER);
    (*centrality) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GRB_TRY(GrB_Matrix_nrows(&n, G->A));

    bool use_all_nodes = (sources == NULL);
    if (!use_all_nodes)
    {
        GRB_TRY(GrB_Vector_nvals(&source_count, sources));
        use_all_nodes = (source_count == 0);
    }

    //--------------------------------------------------------------------------
    // select shortest-path algorithm
    //--------------------------------------------------------------------------

    cc_algo_t algo;
    double emin = -1;

    if (use_floyd_warshall)
    {
        LG_ASSERT_MSG(use_all_nodes, GrB_INVALID_VALUE,
                      "use_floyd_warshall requires sources to be NULL or empty");

        if (!use_weights)
        {
            // Silently downgrade to BFS for unweighted graphs
            // TODO: revisit this? make it more explicit to user?
            algo = CC_BFS;
        }
        else
        {
            // Reject graphs with known or possible negative edges.
            if (G->emin != NULL)
            {
                double fw_emin = 0;
                GRB_TRY(GrB_Scalar_extractElement_FP64(&fw_emin,
                                                       G->emin));
                LG_ASSERT_MSG(fw_emin >= 0, GrB_INVALID_VALUE,
                              "Floyd-Warshall does not support negative edge"
                              " weights");
            }
            algo = CC_FLOYD_WARSHALL;
        }
    }
    else if (!use_weights)
    {
        // Unweighted: BFS is the most efficient choice.
        algo = CC_BFS;
    }
    else
    {
        LG_ASSERT_MSG(G->emin != NULL &&
                          (G->emin_state == LAGraph_VALUE ||
                           G->emin_state == LAGraph_BOUND),
                      LAGRAPH_NOT_CACHED, "G->emin is required");
        GRB_TRY(GrB_Scalar_extractElement_FP64(&emin, G->emin));

        // Use SSSP if no negative edges, Bellman-Ford otherwise
        algo = (emin >= 0 && G->emin_state == LAGraph_VALUE)
                   ? CC_SSSP
                   : CC_BELLMAN_FORD;
    }

    // set up G_AT: graph over the incoming adjacency. Running any shortest-path
    // algorithm from v on G_AT traverses incoming edges, giving distances to v in G.

    bool is_directed = (G->kind != LAGraph_ADJACENCY_UNDIRECTED &&
                        G->is_symmetric_structure != LAGraph_TRUE);

    if (is_directed)
    {
        LG_ASSERT_MSG(G->AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required");
        GRB_TRY(GrB_Matrix_dup(&A_work, G->AT));
        LG_TRY(LAGraph_New(&G_AT, &A_work,
                           LAGraph_ADJACENCY_DIRECTED, msg));

        if (G->emin != NULL)
        {
            GRB_TRY(GrB_Scalar_dup(&G_AT->emin, G->emin));
            G_AT->emin_state = G->emin_state;
        }
    }

    // G_trav->A is the incoming adjacency used by all algorithm paths.
    LAGraph_Graph G_trav = (G_AT != NULL) ? G_AT : G;

    //--------------------------------------------------------------------------
    // build source node list
    //--------------------------------------------------------------------------

    if (!use_all_nodes)
    {
        LG_TRY(LAGraph_Malloc((void **)&source_indices,
                              source_count, sizeof(GrB_Index), msg));
        GRB_TRY(GrB_Vector_extractTuples_UINT64(source_indices, NULL,
                                         &source_count, sources));
        for (GrB_Index k = 0; k < source_count; k++)
        {
            LG_ASSERT(source_indices[k] < n, GrB_INVALID_INDEX);
        }
    }
    else
    {
        source_count = n;
    }

    //--------------------------------------------------------------------------
    // allocate output vector
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Vector_new(&centrality_vector, GrB_FP64, n));
    if (use_all_nodes)
    {
        // Dense output: isolated nodes keep score 0.
        GRB_TRY(GrB_assign(centrality_vector, NULL, NULL,
                           0.0, GrB_ALL, n, NULL));
    }

    //==========================================================================
    // Floyd-Warshall path
    //==========================================================================

    if (algo == CC_FLOYD_WARSHALL)
    {
        GrB_Type D_type;
        GRB_TRY(LAGraph_FW(G_trav->A, &D_APSP, &D_type));

        // Remove self-loops so they do not bias sums or counts.
        GRB_TRY(GrB_select(D_APSP, NULL, NULL,
                           GrB_OFFDIAG, D_APSP, (int64_t)0, NULL));

        // dist_sums[v] = sum of row v
        GRB_TRY(GrB_Vector_new(&dist_sums, GrB_FP64, n));
        GRB_TRY(GrB_reduce(dist_sums, NULL, NULL,
                           GrB_PLUS_MONOID_FP64, D_APSP, NULL));

        // reachable_counts[v] = number of non-zeros in row v of D_APSP.
        // Multiply D_APSP by a dense all-zeros vector with LAGraph_plus_one_fp64:                                                                            
        // TODO: revisit?                                                                            
        GRB_TRY(GrB_Vector_new(&x, GrB_FP64, n));                              
        GRB_TRY(GrB_assign(x, NULL, NULL, (double)0, GrB_ALL, n, NULL));                                                                                      
        GRB_TRY(GrB_Vector_new(&reachable_counts, GrB_FP64, n));                                                                                              
        GRB_TRY(GrB_mxv(reachable_counts, NULL, NULL, LAGraph_plus_one_fp64,
                        D_APSP, x, NULL));                                                                                                                    
        GRB_TRY(GrB_free(&D_APSP));                                                                                                                           
                                                                                                                                                            
        // centrality[v] = R(v) / dist_sums[v].                                                                                                               
        GRB_TRY(GrB_eWiseMult(centrality_vector, NULL, NULL, GrB_DIV_FP64,
                            reachable_counts, dist_sums, NULL));

        (*centrality) = centrality_vector;
        LG_FREE_WORK;
        return GrB_SUCCESS;
    }

    if (algo == CC_SSSP)
    {
        GRB_TRY(GrB_Scalar_new(&Delta, GrB_FP64));
        GRB_TRY(GrB_Scalar_setElement_FP64(Delta, emin > 0 ? emin : 1.0)); // TODO: how to set delta properly?
        GRB_TRY(GrB_Scalar_new(&inf_scalar, GrB_FP64));
        GRB_TRY(GrB_Scalar_setElement_FP64(inf_scalar, (double)INFINITY));
    }

    //==========================================================================
    // per-node shortest-path loop
    //==========================================================================

    for (GrB_Index k = 0; k < source_count; k++)
    {
        GrB_Index node_to_score =
            use_all_nodes ? k : source_indices[k];

        double distance_sum = 0;
        GrB_Index reachable_count = 0;

        //----------------------------------------------------------------------
        // BFS (unweighted shortest paths)
        //----------------------------------------------------------------------

        if (algo == CC_BFS)
        {
            // IMPORTANT: LAGr_BreadthFirstSearch sets *level = NULL without
            // freeing the prior handle.  We must free bfs_level before every
            // call to avoid a memory leak. 
            GRB_TRY(GrB_free(&bfs_level));
            LG_TRY(LAGr_BreadthFirstSearch(&bfs_level, NULL,
                                           G_trav, node_to_score, msg));

            GRB_TRY(GrB_Vector_removeElement(bfs_level, node_to_score));
            GRB_TRY(GrB_Vector_nvals(&reachable_count, bfs_level));
            if (reachable_count > 0)
            {
                GRB_TRY(GrB_reduce(&distance_sum, NULL,
                                   GrB_PLUS_MONOID_FP64, bfs_level, NULL));
            }
        }

        //----------------------------------------------------------------------
        // Delta-stepping SSSP (Dijkstra-like, non-negative weights)
        //----------------------------------------------------------------------

        else if (algo == CC_SSSP)
        {
            GRB_TRY(GrB_free(&distance_vector));
            LG_TRY(LAGr_SingleSourceShortestPath(&distance_vector,
                                                 G_trav,
                                                 node_to_score,
                                                 Delta, msg));

            // Prune unreachable nodes (value == INFINITY).
            GRB_TRY(GrB_select(distance_vector, NULL, NULL,
                               GrB_VALUELT_FP64,
                               distance_vector, inf_scalar, NULL));

            // Exclude self (distance = 0).
            GRB_TRY(GrB_Vector_removeElement(distance_vector,
                                             node_to_score));
            GRB_TRY(GrB_Vector_nvals(&reachable_count, distance_vector));

            if (reachable_count > 0)
            {
                GRB_TRY(GrB_reduce(&distance_sum, NULL,
                                   GrB_PLUS_MONOID_FP64,
                                   distance_vector, NULL));
            }
        }

        //----------------------------------------------------------------------
        // Bellman-Ford 
        //----------------------------------------------------------------------

        else // CC_BELLMAN_FORD
        {
            GRB_TRY(GrB_free(&distance_vector));

            info = LAGraph_BF_basic(&distance_vector,
                                    G_trav->A, node_to_score);
            if (info == GrB_NO_VALUE)
            {
                // Negative-weight cycle reachable from this node;
                // TODO: how to deal with this?
                // For now, set centrality to NaN to indicate an invalid score, and continue.
                GRB_TRY(GrB_Vector_setElement(centrality_vector,
                                              (double)NAN, node_to_score));
                continue;
            }
            GRB_TRY(info);

            // Exclude self (distance = 0).
            GRB_TRY(GrB_Vector_removeElement(distance_vector,
                                             node_to_score));
            GRB_TRY(GrB_Vector_nvals(&reachable_count, distance_vector));

            if (reachable_count > 0)
            {
                GRB_TRY(GrB_reduce(&distance_sum, NULL,
                                   GrB_PLUS_MONOID_FP64,
                                   distance_vector, NULL));
            }
        }

        // closeness(v) = R(v) / sum_{u->v} d(u, v)
        if (reachable_count > 0 && distance_sum > 0)
        {
            double closeness = ((double)reachable_count) / distance_sum;
            GRB_TRY(GrB_Vector_setElement(centrality_vector,
                                          closeness, node_to_score));
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*centrality) = centrality_vector;
    LG_FREE_WORK;
    return GrB_SUCCESS;
}