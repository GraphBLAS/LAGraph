//------------------------------------------------------------------------------
// LAGraph_Leiden.c: community detection using the Leiden algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

#include "GraphBLAS.h"
#include "LG_internal.h"
#include <LAGraphX.h>
#include <LAGraph.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <math.h>

#define LEIDEN_MAX_ITER 20

#ifndef LG_LEIDEN_TIMING
#define LG_LEIDEN_TIMING 1
#endif

#if LG_LEIDEN_TIMING
#define LG_LEIDEN_TIC(t) double t = LAGraph_WallClockTime ( )
#define LG_LEIDEN_ELAPSED(t) (LAGraph_WallClockTime ( ) - (t))
#define LG_LEIDEN_PRINTF(...) printf (__VA_ARGS__)
#define LG_LEIDEN_BURBLE_ON  GRB_TRY (LG_SET_BURBLE (false))
#define LG_LEIDEN_BURBLE_OFF GRB_TRY (LG_SET_BURBLE (false))
#else
#define LG_LEIDEN_TIC(t)
#define LG_LEIDEN_ELAPSED(t) (0.0)
#define LG_LEIDEN_PRINTF(...)
#define LG_LEIDEN_BURBLE_ON
#define LG_LEIDEN_BURBLE_OFF
#endif

#define LG_QUEUE_ENQUEUE(queue, queue_tail, queue_size, value)           \
{                                                                        \
    (queue) [queue_tail] = (value);                                      \
    (queue_tail) = ((queue_tail) + 1) % (queue_size);                    \
}

#define LG_QUEUE_DEQUEUE(queue, queue_head, queue_size, value)           \
{                                                                        \
    (value) = (queue) [queue_head];                                      \
    (queue_head) = ((queue_head) + 1) % (queue_size);                    \
}

//------------------------------------------------------------------------------
// The Leiden algorithm is a modularity-based community detection method that
// guarantees well-connected communities by introducing a Refinement phase
// between the Local-Move and Aggregation phases of the Louvain algorithm.
//
// Reference:
//   Traag, V.A., Waltman, L. & van Eck, N.J. (2019). From Louvain to Leiden:
//   guaranteeing well-connected communities. Scientific Reports 9, 5233.
//   https://doi.org/10.1038/s41598-019-41695-z
//
// Algorithm (one outer iteration = one level of the hierarchy):
//
//   Phase 1 (Local Move): Greedily assign each node to the neighboring
//     community c that maximises the score:
//       score(i->c) = T[c] - k[i] * k_comm[c] / (2m)
//     where T[c] = sum_{j in c, j!=i} A[i,j], k_comm[c] = total degree of
//     community c (excluding i), and m = total edge weight / 2 (constant).
//     This is the standard Louvain/Leiden modularity-gain formula.
//     Initial partition is induced by the Phase-1 communities of the previous
//     level (singletons on the first level).
//     Repeat until no node changes community.
//
//   Phase 2 (Refinement – key Leiden addition): Save Phase-1 communities as
//     "parent" communities.  Restart each node in its own singleton
//     sub-community.  In each local-move step a node may only join a
//     sub-community whose parent equals its own Phase-1 parent.  This ensures
//     every output community is a connected subgraph of the corresponding
//     Phase-1 community.
//
//     NOTE: Traag, Waltman & van Eck (2019) define refinement as a randomized
//     procedure (sampling moves with probability ~ exp(dQ/theta) over moves
//     with dQ >= 0).  This implementation uses the simpler greedy variant: a
//     node is moved to the neighboring sub-community that maximises dQ.  This
//     still satisfies Leiden's connectedness property in practice (every
//     refined sub-community is induced by edges within a Phase-1 community)
//     but does not provide the formal well-connectedness guarantee of the
//     randomized version.  The `seed` parameter is reserved for a future
//     randomized refinement step and is currently unused.
//     Repeat until no node changes sub-community.
//
//   Phase 3 (Aggregation): Build the coarsened graph
//       A_agg = S^T * A_cur * S
//     where S is the n_cur x K_ref membership matrix from Phase 2.  Each
//     column of S is a refined sub-community; each row of A_agg is one
//     super-node.  m is invariant under this operation.  Repeat the outer
//     loop on A_agg until no further coarsening occurs.
//
// Input:  G  – undirected graph (or directed with symmetric structure)
// Output: c_handle – GrB_Vector of INT64 where c[i] = community of node i,
//         labeled 0..K-1.

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free (&c_deg_vec);                                  \
    GrB_free (&c_size_vec);                                 \
    GrB_free (&neighbor_it);                                \
    GrB_free (&it);                                         \
    GrB_free (&x);                                          \
    LAGraph_Free ((void **) &c_deg, NULL) ;                 \
    LAGraph_Free ((void **) &node_deg, NULL) ;              \
    LAGraph_Free ((void **) &community_size, NULL) ;        \
    LAGraph_Free ((void **) &empty_stack, NULL) ;           \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
}

// helper function: phase 1 of leiden. Output already allocated
int LG_Leiden_move_nodes
(
    // output:
    GrB_Matrix C,         // community matrix
    uint64_t *community,  // community array (inside of commmunity matrix)
    uint64_t *nodes_popped_handle,   // queue pops in this phase
    uint64_t *nodes_evaluated_handle, // non-singleton nodes evaluated
    uint64_t *nodes_moved_handle,    // nodes moved to a new community
    // input:
    const GrB_Matrix A,    // adjacency matrix
    const GrB_Vector deg,  // degree vector
    uint64_t *queue,       // queue to use (contains all nodes)
    bool *enqueued,        // nodes in queue (all at the start)
    double m_inv2,       // -1 / (2 * m)
    char* msg
) {
    uint64_t node_id, com_id, n ;
    GrB_Vector x = NULL ;
    GrB_Vector c_deg_vec  = NULL ;
    GrB_Vector c_size_vec = NULL ;
    double *c_deg = NULL ;
    double *node_deg = NULL ;
    uint64_t *community_size = NULL ;
    uint64_t *empty_stack = NULL ;
    uint64_t empty_top = 0 ;
    GxB_Iterator it = NULL, neighbor_it = NULL ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    uint64_t queue_head = 0, queue_tail = n, queue_size = n + 1;

    GRB_TRY (GxB_Iterator_new (&it)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&c_deg_vec, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&c_size_vec, GrB_UINT64, n)) ;
    GRB_TRY (GxB_Iterator_new (&neighbor_it)) ;
    GRB_TRY (GxB_rowIterator_attach (neighbor_it, A, NULL)) ;
    LG_TRY (LAGraph_Malloc ((void **) &empty_stack, n, sizeof (uint64_t), msg)) ;

    // initialize degree and community degree vectors, then unpack to C arrays
    GRB_TRY (GrB_assign (c_deg_vec, NULL, NULL, 0.0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_assign (c_size_vec, NULL, NULL, (uint64_t) 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_vxm (c_deg_vec, NULL, GrB_PLUS_FP64, GxB_PLUS_FIRST_FP64, deg, C, NULL)) ;
    GRB_TRY (GrB_vxm (c_size_vec, NULL, GrB_PLUS_UINT64, GxB_PLUS_PAIR_UINT64, deg, C, NULL)) ;

    GrB_Type node_type = NULL, cdeg_type = NULL, csize_type = NULL ;
    uint64_t node_nheld = 0, cdeg_nheld = 0, csize_nheld = 0 ;
    uint64_t node_size = 0, cdeg_size = 0, csize_size = 0 ;
    int node_handling = 0, cdeg_handling = 0, csize_handling = 0 ;
    GRB_TRY (GxB_Vector_unload (deg, (void **) &node_deg, &node_type,
        &node_nheld, &node_size, &node_handling, NULL)) ;
    GRB_TRY (GxB_Vector_unload (c_deg_vec, (void **) &c_deg, &cdeg_type,
        &cdeg_nheld, &cdeg_size, &cdeg_handling, NULL)) ;
    GRB_TRY (GxB_Vector_unload (c_size_vec, (void **) &community_size, &csize_type,
        &csize_nheld, &csize_size, &csize_handling, NULL)) ;

    for (uint64_t c = 0 ; c < n ; c++)
    {
        if (community_size [c] == 0)
        {
            empty_stack [empty_top++] = c ;
        }
    }

    // TODO: decide if queue should be made inside this function.
    uint64_t nodes_popped = 0, nodes_evaluated = 0, nodes_moved = 0;
    while (queue_head != queue_tail) { // while queue not empty
        nodes_popped++;
        LG_QUEUE_DEQUEUE (queue, queue_head, queue_size, node_id) ;
        enqueued [node_id] = false ;
        uint64_t n_neighbors;
        com_id = community [node_id] ;
        const uint64_t old_size = community_size [com_id] ;
        const double degree_i = node_deg [node_id] ;
        GRB_TRY (GrB_Vector_clear (x));
        // GRB_TRY (GrB_Vector_setElement_FP64 (x, 0.0, node_id)) ;
        // GRB_TRY (GrB_vxm (x, NULL, GrB_PLUS_FP64, GxB_PLUS_SECOND_FP64, x, A, NULL)) ;
        // TODO: deg and c_deg as C_arrays?
        GRB_TRY (GrB_Vector_setElement_FP64 (x, 0.0, node_id)) ;
        GRB_TRY (GrB_Col_extract (x, NULL, GrB_PLUS_FP64, A, GrB_ALL, n, node_id, GrB_DESC_T0)) ;

        GRB_TRY (GrB_Vector_nvals (&n_neighbors, x)) ;
        if (n_neighbors == 1) continue; // skip singletons
        nodes_evaluated++;
        // give gain the number of edges connecting x to community i
        GRB_TRY (GrB_vxm (x, NULL, NULL, GxB_PLUS_FIRST_FP64, x, C, NULL)) ;

        // momentarily remove node degree from the community total
        c_deg [com_id] -= degree_i ;
        community_size [com_id]-- ;

        // Calculate gain in each neighboring community
        uint64_t max_gain_c = com_id ;
        double max_gain_val = -1.0 ;
        double y = m_inv2 * degree_i ;
        GRB_TRY (GxB_Vector_Iterator_attach (it, x, NULL)) ;
        GrB_Info info = GxB_Vector_Iterator_seek (it, 0) ;
        while (info == GrB_SUCCESS) {
            uint64_t c = GxB_Vector_Iterator_getIndex (it) ;
            double gain = GxB_Iterator_get_FP64 (it) ;
            gain += y * c_deg [c] ;
            if (gain > max_gain_val) {
                max_gain_val = gain ;
                max_gain_c = c ;
            }
            info = GxB_Vector_Iterator_next (it) ;
        }

        uint64_t new_com_id = com_id ;
        if (max_gain_val >= 0.0)
        {
            new_com_id = max_gain_c ;
        }
        else if (old_size > 1 && empty_top > 0)
        {
            new_com_id = empty_stack [--empty_top] ;
        }

        community [node_id] = new_com_id ;
        c_deg [new_com_id] += degree_i ;
        community_size [new_com_id]++ ;

        if (new_com_id == com_id) continue; // no change in community
        nodes_moved++ ;
        if (community_size [com_id] == 0)
        {
            empty_stack [empty_top++] = com_id ;
        }

        // push every neighbor in a different community that is not already in
        // the queue to the back of the queue
        info = GxB_rowIterator_seekRow (neighbor_it, node_id) ;
        while (info == GrB_SUCCESS) {
            uint64_t neighbor_id = GxB_rowIterator_getColIndex (neighbor_it) ;
            if (new_com_id != community [neighbor_id]
                && !enqueued [neighbor_id]) {
                LG_QUEUE_ENQUEUE (queue, queue_tail, queue_size, neighbor_id) ;
                enqueued[neighbor_id] = true ;
            }

            info = GxB_rowIterator_nextCol (neighbor_it) ;
        }
    }

    GRB_TRY (GxB_Vector_load (deg, (void **) &node_deg, node_type,
        node_nheld, node_size, node_handling, NULL)) ;

    if (nodes_popped_handle != NULL) (*nodes_popped_handle) = nodes_popped;
    if (nodes_evaluated_handle != NULL) (*nodes_evaluated_handle) = nodes_evaluated;
    if (nodes_moved_handle != NULL) (*nodes_moved_handle) = nodes_moved;
    LG_FREE_WORK ;
    return GrB_SUCCESS ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free (&gain_op);                                    \
    GrB_free (&ctx_s);                                      \
    GrB_free (&leiden_ctx_t);                               \
    GrB_free (&X);                                          \
    GrB_free (&C_t);                                        \
    GrB_free (&it);                                         \
    LAGraph_Free ((void **) &s_deg, NULL) ;                 \
    LAGraph_Free ((void **) &is_rep, NULL) ;                \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
}

typedef struct {
    double *s_deg;
    double km_inv;
} leiden_ctx ;

void LG_Leiden_gain
(
    double *z,
    const double *x,
    GrB_Index i,
    GrB_Index j,
    const leiden_ctx *thunk
) {
    *z = *x + thunk->km_inv * thunk->s_deg [j];
}
// helper function: phase 2 of leiden. Output already allocated
int LG_Leiden_refinement
(
    // output:
    GrB_Matrix S,         // sub-community matrix
    uint64_t *sub_com,    // sub-community array (inside of S matrix)
    uint64_t *nodes_evaluated_handle, // non-singleton nodes evaluated
    // input:
    const GrB_Matrix C,         // community matrix
    const uint64_t *community,  // community array (inside of C matrix)
    const GrB_Matrix A,         // adjacency matrix
    const GrB_Vector deg,       // degree vector
    uint64_t *queue,            // queue to use (contains all nodes)
    double m_inv2,              // -1 / (2 * m)
    char* msg
) {
    uint64_t node_id, com_id, n_deg, n;
    GrB_Matrix X = NULL;
    double *s_deg = NULL;
    bool *is_rep = NULL;
    GrB_Matrix C_t = NULL ; // transpose of C
    GxB_Iterator it = NULL ;
    GrB_IndexUnaryOp gain_op = NULL;
    GrB_Type leiden_ctx_t = NULL;
    GrB_Scalar ctx_s = NULL;

    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    uint64_t queue_head = 0, queue_tail = n, queue_size = n + 1;

    // TODO: make C_t compact, use an different queue per each community.
    // This makes it easier to do the per-comunity masking.

    // TODO: decide if queue should be made inside this function.
    GRB_TRY (GrB_Matrix_new (&C_t, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_transpose (C_t, NULL, NULL, C, NULL)) ;

    uint64_t n_communities;

    LG_TRY (LAGraph_Malloc ((void **) &is_rep, n, sizeof (bool), msg)) ;
    memset (is_rep, 0, n) ;
    LG_TRY (LAGraph_Malloc ((void **) &s_deg, n, sizeof (uint64_t), msg)) ;
    GRB_TRY (GrB_Vector_extractTuples_FP64 (NULL, s_deg, &n, deg)) ;
    GRB_TRY (GrB_Type_new (&leiden_ctx_t, sizeof(leiden_ctx))) ;
    GRB_TRY (GrB_Matrix_new (&X, GrB_FP64, n, n)) ;
    // Work around SuiteSparse:GraphBLAS bug (67) present in v10.3.0 and
    // fixed in v10.3.1: incorrect JIT kernel for R=masker(C,M,Z) when R is
    // hypersparse. Apply only to affected versions by disabling hypersparsity
    // for X on this path.
    #if LAGRAPH_SUITESPARSE && (GxB_IMPLEMENTATION < GxB_VERSION (10,3,1))
    GRB_TRY (GxB_Matrix_Option_set_FP64 (X, GxB_HYPER_SWITCH, 0.0)) ;
    #endif
    GRB_TRY (GrB_set (X, GxB_SPARSE, GxB_SPARSITY_CONTROL)) ;
    GRB_TRY (GrB_Scalar_new (&ctx_s, leiden_ctx_t)) ;
    GRB_TRY (GrB_IndexUnaryOp_new (
        &gain_op, (GxB_index_unary_function) LG_Leiden_gain,
        GrB_FP64, GrB_FP64, leiden_ctx_t)) ;
    GRB_TRY (GxB_Iterator_new(&it)) ;
    leiden_ctx ctx = {.s_deg = s_deg, .km_inv = 0.0} ;

    uint64_t nodes_evaluated = 0;
    while (queue_head != queue_tail) { // while queue not empty
        double count;
        LG_QUEUE_DEQUEUE (queue, queue_head, queue_size, node_id) ;

        if (is_rep [node_id]) {
            continue;
        }
        nodes_evaluated++;

        double com_deg, n_deg;
        com_id = community [node_id] ;
        uint64_t sub_com_id = sub_com [node_id] ;

        // TODO: this is probably much easier with extract and if we select
        // nodes via their communities.
        GRB_TRY (GrB_Matrix_clear (X));
        GRB_TRY (GrB_Matrix_setElement_FP64 (X, 0.0, com_id, node_id));
        GRB_TRY (GrB_mxm (X, C_t, NULL, GxB_PLUS_SECOND_FP64, X, A, GrB_DESC_S)) ;
        // give X the number of edges connecting x to community i
        GRB_TRY (GrB_mxm (X, NULL, NULL, GxB_PLUS_FIRST_FP64, X, S, NULL)) ;

        // momentarily remove node degree from the community total
        GRB_TRY (GrB_Vector_extractElement_FP64 (&n_deg, deg, node_id));

        s_deg [sub_com_id] -= n_deg;

        ctx.km_inv = m_inv2 * n_deg;
        GRB_TRY (GrB_Scalar_setElement_UDT (ctx_s, &ctx)) ;
        // Calculate gain in each neighboring community.
        GRB_TRY (GrB_apply (X, NULL, NULL, gain_op, X, ctx_s, NULL)) ;

        // FUTURE: This part will have to change if we want parrallel clusters
        // to work
        uint64_t new_c = sub_com_id, row = 0;
        // Select the best gain candidate (temporary greedy refinement).
        GRB_TRY (GxB_Matrix_Iterator_attach(it, X, NULL)) ;
        GrB_Info info = GxB_Matrix_Iterator_seek (it, 0) ;
        double best_gain = 0.0 ;
        while (info == GrB_SUCCESS) {
            double gain = GxB_Iterator_get_FP64 (it) ;
            if (gain > best_gain) {
                best_gain = gain ;
                GxB_Matrix_Iterator_getIndex(it, &row, &new_c) ;
            }

            info = GxB_Matrix_Iterator_next (it);
        }

        is_rep [new_c] = true;
        sub_com [node_id] = new_c;
        s_deg [new_c] += n_deg;
    }
    if (nodes_evaluated_handle != NULL) (*nodes_evaluated_handle) = nodes_evaluated;
    LG_FREE_WORK ;
    return GrB_SUCCESS ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free (&row_parent_count);                           \
    GrB_free (&row_split_count);                            \
    GrB_free (&all_ones);                                   \
    GrB_free (&C_tS_in_C_t);                                \
    GrB_free (&C_tS);                                       \
    GrB_free (&C_t);                                        \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
}

// helper function: phase 2 metrics (split counts + parent containment check)
int LG_Leiden_refinement_stats
(
    // output:
    uint64_t *split_total_handle, // total extra splits beyond 1 per parent
    double *split_avg_handle,     // average extra splits per parent
    bool *within_parent_handle,   // true iff C_tS is a submatrix of C_t
    // input:
    const GrB_Matrix C,           // parent-community matrix
    const GrB_Matrix S,           // refined sub-community matrix
    char *msg
)
{
    (void) msg;
    GrB_Index n;
    GrB_Matrix C_t = NULL, C_tS = NULL, C_tS_in_C_t = NULL;
    GrB_Vector all_ones = NULL, row_split_count = NULL, row_parent_count = NULL;

    GRB_TRY (GrB_Matrix_nrows (&n, C)) ;
    GRB_TRY (GrB_Matrix_new (&C_t, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_transpose (C_t, NULL, NULL, C, NULL)) ;
    GRB_TRY (GrB_assign (
        C_t, C_t, NULL, (bool) true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;

    GRB_TRY (GrB_Matrix_new (&C_tS, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_mxm (C_tS, NULL, NULL, LAGraph_any_one_bool, C_t, S, NULL)) ;
    GRB_TRY (GrB_assign (
        C_tS, C_tS, NULL, (bool) true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;

    GRB_TRY (GrB_Vector_new (&all_ones, GrB_BOOL, n)) ;
    GRB_TRY (GrB_assign (all_ones, NULL, NULL, (bool) true, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_Vector_new (&row_split_count, GrB_INT32, n)) ;
    GRB_TRY (GrB_Vector_new (&row_parent_count, GrB_INT32, n)) ;
    GRB_TRY (GrB_mxv (
        row_split_count, NULL, NULL, LAGraph_plus_one_int32, C_tS, all_ones, NULL)) ;
    GRB_TRY (GrB_mxv (
        row_parent_count, NULL, NULL, LAGraph_plus_one_int32, C_t, all_ones, NULL)) ;

    int64_t split_total_raw = 0;
    GRB_TRY (GrB_reduce (
        &split_total_raw, NULL, GrB_PLUS_MONOID_INT64, row_split_count, NULL)) ;
    GrB_Index n_parent = 0;
    GRB_TRY (GrB_Vector_nvals (&n_parent, row_parent_count)) ;
    int64_t split_total = split_total_raw - (int64_t) n_parent;
    if (split_total < 0) split_total = 0;

    if (split_total_handle != NULL) (*split_total_handle) = (uint64_t) split_total;
    if (split_avg_handle != NULL) {
        (*split_avg_handle) =
            (n_parent == 0) ? 0.0 : ((double) split_total / (double) n_parent);
    }

    // Verify that every (parent, sub-community) incidence remains in C_t.
    GRB_TRY (GrB_Matrix_new (&C_tS_in_C_t, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_eWiseMult (
        C_tS_in_C_t, NULL, NULL, GrB_LAND, C_tS, C_t, NULL)) ;
    GrB_Index nvals_C_tS = 0, nvals_intersection = 0;
    GRB_TRY (GrB_Matrix_nvals (&nvals_C_tS, C_tS)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals_intersection, C_tS_in_C_t)) ;
    if (within_parent_handle != NULL) {
        (*within_parent_handle) = (nvals_C_tS == nvals_intersection);
    }

    LG_FREE_WORK ;
    return GrB_SUCCESS ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                         \
{                                                            \
    GrB_free (&desc) ;                                       \
    GrB_free (&x) ;                                          \
    GrB_free (&s_list) ;                                     \
    GrB_free (&S_squished) ;                                 \
    GrB_free (&C_squished) ;                                 \
    GrB_free (&A_squished) ;                                 \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                    \
{                                                      \
    GrB_free (S) ;                                     \
    GrB_free (C) ;                                     \
    if (free_A) GrB_free (A) ;                         \
    LG_FREE_WORK ;                                     \
}

// helper function: phase 3 of leiden. input freed and output allocated
int LG_Leiden_aggregate
(
    // input / output:
    GrB_Matrix *S,         // sub-communities to agregate
    GrB_Matrix *C,         // communities to preserve
    GrB_Matrix *A,         // adjacency matrix
    bool free_A,           // can A be freed by this function?
    char* msg
) {
    uint64_t node_id, com_id, n_deg, queue_head = 0, queue_tail = 0;
    GrB_Vector s_list = NULL, x = NULL;
    // TODO:rename S_new, maybe make it pure output instead of i/o
    GrB_Matrix S_squished = NULL, S_new = NULL ;
    GrB_Matrix C_squished = NULL, C_new = NULL ;
    GrB_Matrix A_squished = NULL, A_new = NULL ;

    GrB_Descriptor desc = NULL ;
    uint64_t n, n_new ;
    GRB_TRY (GrB_Matrix_nrows (&n, *A)) ;

    GRB_TRY (GrB_Descriptor_new (&desc)) ;
    GRB_TRY (GrB_set (desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)) ;
    GRB_TRY (GrB_set (desc, GxB_USE_INDICES, GxB_COLINDEX_LIST)) ;

    GRB_TRY (GrB_Vector_new (&s_list, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_BOOL, n)) ;
    GRB_TRY (GrB_assign (x, NULL, NULL, (bool) 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_vxm (s_list, NULL, NULL, GxB_ANY_PAIR_BOOL, x, *S, NULL)) ;

    GRB_TRY (GrB_Vector_nvals (&n_new, s_list)) ;

    GRB_TRY (GrB_Matrix_new (&S_new, GrB_BOOL, n_new, n)) ;
    GRB_TRY (GrB_Matrix_new (&C_new, GrB_BOOL, n_new, n_new)) ;
    GRB_TRY (GrB_Matrix_new (&A_new, GrB_FP64, n_new, n_new)) ;
    GrB_set (A_new, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT) ;

    GRB_TRY (GrB_Matrix_new (&S_squished, GrB_BOOL, n, n_new)) ;
    GRB_TRY (GrB_Matrix_new (&C_squished, GrB_BOOL, n, n_new)) ;
    GRB_TRY (GrB_Matrix_new (&A_squished, GrB_FP64, n, n_new)) ;

    GRB_TRY (GxB_Matrix_extract_Vector (S_squished, NULL, NULL, *S, NULL, s_list, desc)) ;
    GRB_TRY (GrB_free (S)) ;

    GRB_TRY (GrB_mxm (C_squished, NULL, NULL, GxB_ANY_PAIR_BOOL, *C, S_squished, NULL)) ;
    GRB_TRY (GrB_free (C)) ;

    GRB_TRY (GrB_mxm (A_squished, NULL, NULL, GxB_PLUS_FIRST_FP64, *A, S_squished, NULL)) ;
    if (free_A)
        GRB_TRY (GrB_free (A)) ;

    GRB_TRY (GrB_transpose (S_new, NULL, NULL, S_squished, NULL)) ;
    GRB_TRY (GrB_free (&S_squished)) ;

    GRB_TRY (GrB_mxm (C_new, NULL, NULL, GxB_ANY_PAIR_BOOL, S_new, C_squished, NULL)) ;
    GRB_TRY (GrB_free (&C_squished)) ;

    GRB_TRY (GrB_mxm (A_new, NULL, NULL, GxB_PLUS_SECOND_FP64, S_new, A_squished, NULL)) ;
    GRB_TRY (GrB_free (&A_squished)) ;

    *S = S_new ;
    *C = C_new ;
    *A = A_new ;

    LG_FREE_WORK;
    return GrB_SUCCESS ;
}
#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free (&neighbor_it);                                \
    GrB_free (&cont);                                       \
    GrB_free (&C_to_orig);                                  \
    GrB_free (&C);                                          \
    GrB_free (&S);                                          \
    if (free_A) GrB_free (&A) ;                             \
    GrB_free (&deg);                                        \
    GrB_free (&x);                                          \
    GrB_free (&node_map);                                   \
    GrB_free (&new_node_map);                               \
    LAGraph_Free ((void **) &queue, NULL);                  \
    LAGraph_Free ((void **) &enqueued, NULL);               \
    LAGraph_Free ((void **) &community, NULL);              \
    LAGraph_Free ((void **) &sub_com, NULL);                \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
}

int LAGraph_Leiden
(
    // output:
    GrB_Vector *c_handle,   // c[i] = community label (0..K-1) for node i
    // input:
    const LAGraph_Graph G,  // input graph (must be symmetric, no self-loops,
                            // have positive numerical weights)
    uint64_t seed,          // random seed
    char *msg
)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10
    uint64_t *queue = NULL; // circular buffer with next nodes to look up
    bool *enqueued = NULL; // true if node is in the queue
    uint64_t *community = NULL ; // array: i in community c[i]
    GrB_Matrix C = NULL ; // contains readonly pointer to community array
    GxB_Container cont = NULL;

    uint64_t *sub_com = NULL ; // array: i in sub-community sub_com[i]
    GrB_Matrix S = NULL ; // contains readonly pointer to sub-community array

    GrB_Vector deg = NULL ; // degree of each node
    GrB_Matrix A = NULL ; // G->A
    GrB_Vector x = NULL; // full boolean vector
    GrB_Matrix node_map = NULL ; // maps aggregate node to nodes it contains
    GrB_Matrix new_node_map = NULL ; // maps aggregate node to nodes it contains
    GrB_Matrix C_to_orig = NULL ; // maps final C communities to original nodes

    uint64_t n, original_n;

    GxB_Iterator neighbor_it = NULL;

    bool free_A = false ;
    A = G->A;
    GRB_TRY (GrB_Matrix_nrows (&n, A));
    original_n = n;  // Store original number of nodes

    // Input checking
    LG_ASSERT_MSG (
        G->is_symmetric_structure == LAGraph_TRUE,
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->emin_state != LAGraph_STATE_UNKNOWN, LAGRAPH_NOT_CACHED,
                   "G->emin must be defined") ;
    // cast to FP64 if needed
    double min_val = 0.0;
    GRB_TRY (GrB_Scalar_extractElement_FP64 (&min_val, G->emin)) ;
    LG_ASSERT_MSG  (min_val >= 0.0, GrB_INVALID_VALUE,
                   "G->emin must be non-negative") ;
    // TODO: check numeric, and should we allow 0 weight edges?

    // initialize queue
    LG_TRY (LAGraph_Malloc ((void **) &enqueued, n, sizeof(bool), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &queue, n + 1, sizeof (uint64_t), msg)) ;

    // Initialize GraphBLAS data structures
    GRB_TRY (GxB_Iterator_new(&neighbor_it));
    GRB_TRY (GxB_Container_new (&cont)) ;
    GRB_TRY (GrB_Vector_new(&deg, GrB_FP64, n));

    GRB_TRY (GrB_Vector_new(&x, GrB_BOOL, n));
    GRB_TRY (GrB_assign (x, NULL, NULL, (bool) 0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_Matrix_diag (&C, x, 0)) ;
    GRB_TRY (GrB_Matrix_dup (&S, C)) ;

    // Initialize full degree vector from graph
    GRB_TRY (GrB_assign (deg, NULL, NULL, 0.0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_reduce (deg, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, A, NULL));
    // don't use out degree unless A is bool adj matrix
    // GRB_TRY (GrB_assign (deg, NULL, GrB_PLUS_FP64, G->out_degree, GrB_ALL, n, NULL));


    double m = 0;
    GRB_TRY (GrB_reduce (&m, NULL, GrB_PLUS_MONOID_FP64, deg, NULL)) ;
    LG_ASSERT_MSG (isfinite(m), GrB_INVALID_VALUE,
                   "Matrix must reduce to a finite value") ;
    double m_inv2 = -1.0 / (m);

    srand(seed);
    srand48 ((long int) seed) ;

    // Outer loop: repeat phases until convergence
    double prev_modularity = -1;
    double modularity_epsilon = 1e-6;

    double total_move_time = 0, total_refine_time = 0, total_agg_time = 0;
    double total_modularity_time = 0, total_unpack_time = 0, total_queue_time = 0;
    uint64_t total_move_evaluated = 0, total_move_popped = 0, total_move_moved = 0;
    uint64_t total_refine_evaluated = 0;
    uint64_t total_refine_splits = 0, total_refine_calls = 0;
    double total_refine_split_avg = 0.0;

    for (int count = 0; count < LEIDEN_MAX_ITER; count++) {
        double current_modularity;
        LG_LEIDEN_PRINTF ("[Leiden timing] iteration %d, coarse nodes=%llu\n",
            count, (unsigned long long) n) ;
        // Initialize queue in random order (seed-based)
        LG_LEIDEN_TIC (t_queue_phase1);
        uint64_t queue_len = n;
        // TODO: store values and memcopy in a different array
        GRB_TRY (GrB_Vector_extractTuples_FP64 (queue, NULL, &queue_len, deg)) ;

        // TODO: shuffle with LAGraph Random instead.
        if (queue_len > 1) {
            for (uint64_t i = queue_len - 1; i > 0; i--) {
                uint64_t j = rand() % (i + 1);
                uint64_t temp = queue[i];
                queue[i] = queue[j];
                queue[j] = temp;
            }
        }
        memset (enqueued, true, n) ;
        double queue_phase1_time = LG_LEIDEN_ELAPSED (t_queue_phase1);
        total_queue_time += queue_phase1_time;

        // get i vector from C and S to use as community and sub_com arrays
        // keep inside C and S as read-only
        LG_LEIDEN_TIC (t_unpack);

        {
            GrB_Type type;
            uint64_t n_community, X_memsize;
            int handling;

            GRB_TRY (GrB_set (C, GxB_SPARSE, GxB_SPARSITY_CONTROL)) ;
            GRB_TRY (GrB_set (C, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
            GRB_TRY (GrB_set (C, 64, GxB_COLINDEX_INTEGER_HINT)) ;
            GrB_wait (C, GrB_MATERIALIZE) ;
            GRB_TRY (GxB_unload_Matrix_into_Container(C, cont, NULL)) ;
            GRB_TRY (GxB_Vector_unload (cont->i, (void **) &community, &type, &n_community,
                &X_memsize, &handling, NULL)) ;
            GRB_TRY (GxB_Vector_load (cont->i, (void **) &community, type, n_community,
                X_memsize, GxB_IS_READONLY, NULL)) ;
            GRB_TRY (GxB_load_Matrix_from_Container(C, cont, NULL)) ;

            GRB_TRY (GrB_set (S, GxB_SPARSE, GxB_SPARSITY_CONTROL)) ;
            GRB_TRY (GrB_set (S, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
            GRB_TRY (GrB_set (S, 64, GxB_COLINDEX_INTEGER_HINT)) ;
            GrB_wait (S, GrB_MATERIALIZE) ;
            GRB_TRY (GxB_unload_Matrix_into_Container (S, cont, NULL)) ;
            GRB_TRY (GxB_Vector_unload (cont->i, (void **) &sub_com, &type, &n_community,
                &X_memsize, &handling, NULL)) ;
            GRB_TRY (GxB_Vector_load (cont->i, (void **) &sub_com, type, n_community,
                X_memsize, GxB_IS_READONLY, NULL)) ;
            GRB_TRY (GxB_load_Matrix_from_Container(S, cont, NULL)) ;
        }

        double unpack_time = LG_LEIDEN_ELAPSED (t_unpack);
        total_unpack_time += unpack_time;

        //----------------------------------------------------------------------
        // Phase 1: move nodes
        //----------------------------------------------------------------------
        uint64_t move_popped = 0, move_evaluated = 0, move_moved = 0;
        LG_LEIDEN_TIC (t_move);
        // LG_LEIDEN_BURBLE_ON;
        LG_TRY (LG_Leiden_move_nodes (C, community, &move_popped, &move_evaluated,
            &move_moved, A, deg, queue, enqueued, m_inv2, msg)) ;
        LG_LEIDEN_BURBLE_OFF;
        double move_time = LG_LEIDEN_ELAPSED (t_move);
        total_move_time += move_time;
        total_move_popped += move_popped;
        total_move_evaluated += move_evaluated;
        total_move_moved += move_moved;
        LG_LEIDEN_PRINTF (
            "  move_nodes: t=%.6fs, queue_pops=%llu, evaluated=%llu, moved=%llu, amortized=%.3fus/evaluated\n",
            move_time,
            (unsigned long long) move_popped,
            (unsigned long long) move_evaluated,
            (unsigned long long) move_moved,
            (move_evaluated == 0) ? 0.0 : (1e6 * move_time / (double) move_evaluated)) ;

    #ifndef NDEBUG
        LG_TRY (LAGr_AdjModularity (&current_modularity, 1.0, A, C, msg));
        ASSERT (current_modularity >= prev_modularity - modularity_epsilon) ;
        LG_LEIDEN_PRINTF ( "  modularity %f\n", current_modularity) ;
    #endif

        //----------------------------------------------------------------------
        // Prep phase 2
        //----------------------------------------------------------------------
        // Initialize queue
        LG_LEIDEN_TIC (t_queue_phase2);
        queue_len = n;
        GRB_TRY (GrB_Vector_extractTuples_FP64 (queue, NULL, &queue_len, deg)) ;

        // TODO: shuffle with LAGraph Random instead.
        if (queue_len > 1) {
            for (uint64_t i = queue_len - 1; i > 0; i--) {
                uint64_t j = rand() % (i + 1);
                uint64_t temp = queue[i];
                queue[i] = queue[j];
                queue[j] = temp;
            }
        }
        double queue_phase2_time = LG_LEIDEN_ELAPSED (t_queue_phase2);
        total_queue_time += queue_phase2_time;

        //----------------------------------------------------------------------
        // Phase 2: refine
        //----------------------------------------------------------------------
        // Save Phase 1 communities as parent communities
        uint64_t refine_evaluated = 0;
        uint64_t refine_split_total = 0;
        double refine_split_avg = 0.0;
        bool refine_within_parent = false;
        LG_LEIDEN_TIC (t_refine);
        LG_LEIDEN_BURBLE_ON;
        // memcpy (sub_com, community, n * sizeof (uint64_t)) ;
        LG_TRY (LG_Leiden_refinement (S, sub_com, &refine_evaluated, C, community,
             A, deg, queue, m_inv2, msg)) ;
        LG_TRY (LG_Leiden_refinement_stats (
            &refine_split_total, &refine_split_avg, &refine_within_parent,
            C, S, msg)) ;
        LG_LEIDEN_BURBLE_OFF;
        double refine_time = LG_LEIDEN_ELAPSED (t_refine);
        total_refine_time += refine_time;
        total_refine_evaluated += refine_evaluated;
        total_refine_splits += refine_split_total;
        total_refine_split_avg += refine_split_avg;
        total_refine_calls++;
        LG_LEIDEN_PRINTF (
            "  refinement: t=%.6fs, evaluated=%llu, splits_total=%llu, splits_avg=%.6g, within_parent=%s, amortized=%.3fus/evaluated\n",
            refine_time,
            (unsigned long long) refine_evaluated,
            (unsigned long long) refine_split_total,
            refine_split_avg,
            refine_within_parent ? "true" : "false",
            (refine_evaluated == 0) ? 0.0 : (1e6 * refine_time / (double) refine_evaluated)) ;

        //----------------------------------------------------------------------
        // Phase 3: aggregate
        //----------------------------------------------------------------------

        // Compute modularity to check convergence
        LG_LEIDEN_TIC (t_modularity);
        // LG_LEIDEN_BURBLE_ON;
        LG_TRY (LAGr_AdjModularity (&current_modularity, 1.0, A, C, msg));
        LG_LEIDEN_BURBLE_OFF;
        double modularity_time = LG_LEIDEN_ELAPSED (t_modularity);
        total_modularity_time += modularity_time;
        LG_LEIDEN_PRINTF ("  modularity: Q=%.12g, t=%.6fs\n",
            current_modularity, modularity_time) ;

        ASSERT (current_modularity >= prev_modularity - modularity_epsilon) ;
        if (fabs(current_modularity - prev_modularity) < modularity_epsilon) {
            break;  // converged
        }
        prev_modularity = current_modularity;

        LG_LEIDEN_TIC (t_agg);
        // LG_LEIDEN_BURBLE_ON;
        LG_TRY (LG_Leiden_aggregate (&S, &C, &A, free_A, msg)) ;
        LG_LEIDEN_BURBLE_OFF;
        double agg_time = LG_LEIDEN_ELAPSED (t_agg);
        total_agg_time += agg_time;
        LG_LEIDEN_PRINTF ("  aggregate: t=%.6fs\n", agg_time) ;
        // free internal pointers that were not freed
        LG_TRY (LAGraph_Free ((void **) &community, msg)) ;
        LG_TRY (LAGraph_Free ((void **) &sub_com, msg)) ;
        GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

        free_A = true;
        if (node_map == NULL) {
            node_map = S;
            S = NULL;
        } else {
            GRB_TRY (GrB_Matrix_new (&new_node_map, GrB_BOOL, n, original_n)) ;
            GRB_TRY (GrB_mxm (new_node_map, NULL, NULL, GxB_ANY_PAIR_BOOL, S, node_map, NULL)) ;

            GRB_TRY (GrB_free (&node_map)) ;
            node_map = new_node_map; new_node_map = NULL;
            GRB_TRY (GrB_free (&S)) ;
        }

        // Re-initialize degree from coarsened graph
        GRB_TRY (GrB_free (&deg)) ;
        GRB_TRY (GrB_Vector_new (&deg, GrB_FP64, n)) ;
        GRB_TRY (GrB_assign(deg, NULL, NULL, 0.0, GrB_ALL, n, NULL));
        GRB_TRY (GrB_reduce(deg, NULL, NULL, GrB_PLUS_MONOID_FP64, A, NULL));
        GRB_TRY (GrB_Vector_resize (x, n)) ;
        GRB_TRY (GrB_Matrix_diag (&S, x, 0)) ;
    #ifndef NDEBUG
        LG_TRY (LAGr_AdjModularity (&current_modularity, 1.0, A, C, msg));
        ASSERT (isfinite (current_modularity)) ;
    #endif
    } // end outer while loop

    LG_LEIDEN_PRINTF (
        "[Leiden timing] totals: queue=%.6fs, unpack=%.6fs, move=%.6fs, refine=%.6fs, modularity=%.6fs, aggregate=%.6fs\n",
        total_queue_time, total_unpack_time, total_move_time, total_refine_time,
        total_modularity_time, total_agg_time) ;
    LG_LEIDEN_PRINTF (
        "[Leiden timing] move totals: queue_pops=%llu, evaluated=%llu, moved=%llu, amortized=%.3fus/evaluated\n",
        (unsigned long long) total_move_popped,
        (unsigned long long) total_move_evaluated,
        (unsigned long long) total_move_moved,
        (total_move_evaluated == 0) ? 0.0 :
            (1e6 * total_move_time / (double) total_move_evaluated)) ;
    LG_LEIDEN_PRINTF (
        "[Leiden timing] refinement totals: evaluated=%llu, splits_total=%llu, splits_avg=%.6g, amortized=%.3fus/evaluated\n",
        (unsigned long long) total_refine_evaluated,
        (unsigned long long) total_refine_splits,
        (total_refine_calls == 0) ? 0.0 :
            (total_refine_split_avg / (double) total_refine_calls),
        (total_refine_evaluated == 0) ? 0.0 :
            (1e6 * total_refine_time / (double) total_refine_evaluated)) ;


    // Map communities in C back to original nodes
    GRB_TRY (GrB_Vector_new(c_handle, GrB_INT64, original_n));
    if (node_map == NULL) {
        ASSERT (n == original_n) ;
        GRB_TRY (GxB_Vector_load (*c_handle, (void **) &community, GrB_INT64, n,
            sizeof (int64_t) * n, GrB_DEFAULT, NULL));
    } else {
        GRB_TRY (GrB_Matrix_new (&C_to_orig, GrB_BOOL, n, original_n)) ;
        GRB_TRY (GrB_mxm (
            C_to_orig, NULL, NULL, GxB_ANY_PAIR_BOOL, C, node_map, GrB_DESC_T0)) ;
        GRB_TRY (GrB_Vector_resize (x, n)) ;
        GRB_TRY (GrB_assign (
            *c_handle, NULL, NULL, (int64_t) 0, GrB_ALL, original_n, NULL)) ;
        GRB_TRY (GrB_vxm (*c_handle, NULL, GrB_PLUS_INT64,
            GxB_PLUS_FIRSTJ_INT64, x, C_to_orig, NULL)) ;
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
#else
    return GrB_NOT_IMPLEMENTED;
#endif
}
