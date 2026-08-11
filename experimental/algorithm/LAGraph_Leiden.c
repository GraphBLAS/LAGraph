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
    GrB_free (&it);                                         \
    GrB_free (&c_deg);                                      \
    GrB_free (&x);                                          \
    GrB_free (&gain);                                       \
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
    // input:
    const GrB_Matrix A,    // adjacency matrix
    const GrB_Vector deg,  // degree vector
    uint64_t *queue,       // queue to use (contains all nodes)
    bool *enqueued,        // nodes in queue (all at the start)
    double m_inv2,       // -1 / (2 * m)
    char* msg
) {
    uint64_t node_id, com_id, n_deg, n;
    GrB_Vector x = NULL, gain = NULL;
    GrB_Vector c_deg = NULL ;
    GxB_Iterator it = NULL;
    GrB_Matrix_nrows (&n, A) ;
    uint64_t queue_head = 0, queue_tail = n, queue_size = n + 1;

    GRB_TRY (GxB_Iterator_new (&it)) ;
    GRB_TRY (GrB_Vector_new (&c_deg, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&gain, GrB_FP64, n)) ;

    GRB_TRY (GrB_mxv(c_deg, NULL, NULL, GxB_PLUS_SECOND_FP64, C, deg, NULL));
    // TODO: decide if queue should be made inside this function.
    while (queue_head != queue_tail) { // while queue not empty
        node_id = queue[queue_head];
        queue_head = (queue_head + 1) % queue_size;
        enqueued [node_id] = false;
        double com_deg, n_deg;
        uint64_t n_neighbors;
        com_id = community [node_id] ;
        GRB_TRY (GrB_Vector_clear (x));
        GRB_TRY (GrB_Vector_setElement_FP64 (x, 0.0, node_id));
        GRB_TRY (GrB_vxm (x, NULL, GrB_PLUS_FP64, GxB_PLUS_SECOND_FP64, x, A, NULL)) ;
        GRB_TRY (GrB_Vector_nvals (&n_neighbors, x)) ;
        if (n_neighbors == 1) continue; // skip singletons
        // give gain the number of edges connecting x to community i
        GRB_TRY (GrB_vxm (gain, NULL, NULL, GxB_PLUS_FIRST_FP64, x, C, NULL)) ;

        // momentarily remove node degree from the community total
        GRB_TRY (GrB_Vector_extractElement_FP64 (&n_deg, deg, node_id));

        // c_deg [com_id] -= ndeg
        GRB_TRY (GrB_assign (c_deg, NULL, GrB_MINUS_FP64, n_deg, &com_id, 1, NULL));

        // Calculate gain in each neighboring community
        GRB_TRY (GrB_apply (gain, gain, GrB_PLUS_FP64, GrB_TIMES_FP64, c_deg, m_inv2 * n_deg, GrB_DESC_S));
        uint64_t max_gain_c;
        double max_gain_val = -1e100;
        GRB_TRY (GxB_Vector_Iterator_attach (it, gain, NULL));
        GrB_Info info = GxB_Vector_Iterator_seek (it, 0);
        while (info == GrB_SUCCESS) {
            uint64_t c = GxB_Vector_Iterator_getIndex (it);
            double g = GxB_Iterator_get_FP64 (it);
            if (g > max_gain_val) {
                max_gain_val = g;
                max_gain_c = c;
            }
            info = GxB_Vector_Iterator_next (it);
        }

        community [node_id] = max_gain_c;
        // c_deg [node_id] -= ndeg
        GRB_TRY (GrB_assign (c_deg, NULL, GrB_PLUS_FP64, n_deg, &max_gain_c, 1, NULL));

        if (max_gain_c == com_id) continue; // no change in community

        // push every neighbor in a different community that is not already in
        // the queue to the back of the queue
        GRB_TRY (GxB_Vector_Iterator_attach (it, x, NULL));
        info = GxB_Vector_Iterator_seek (it, 0);
        while (info == GrB_SUCCESS) {
            uint64_t neighbor_id = GxB_Vector_Iterator_getIndex (it);
            if (max_gain_c != community [neighbor_id]
                && !enqueued [neighbor_id]) {
                queue[queue_tail] = neighbor_id;
                queue_tail = (queue_tail + 1) % queue_size;
                enqueued[neighbor_id] = true;
            }

            info = GxB_Vector_Iterator_next (it);
        }
    }
    LG_FREE_WORK ;
    return GrB_SUCCESS ;
}

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free (&s_deg);                                      \
    GrB_free (&x_count);                                    \
    GrB_free (&k);                                          \
    GrB_free (&X);                                          \
    GrB_free (&C_t);                                        \
    GrB_free (&it);                                         \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
}

// helper function: phase 2 of leiden. Output already allocated
int LG_Leiden_refinement
(
    // output:
    GrB_Matrix S,         // sub-community matrix
    uint64_t *sub_com,    // sub-community array (inside of S matrix)
    // input:
    const GrB_Matrix C,         // community matrix
    const uint64_t *community,  // community array (inside of C matrix)
    const GrB_Vector c_deg,     // degree vector
    const GrB_Matrix A,         // adjacency matrix
    const GrB_Vector deg,       // degree vector
    uint64_t *queue,            // queue to use (contains all nodes)
    double m_inv2,              // -1 / (2 * m)
    char* msg
) {
    uint64_t node_id, com_id, n_deg, n;
    GrB_Matrix X = NULL;
    GrB_Vector s_deg = NULL;
    GrB_Vector k = NULL;
    GrB_Matrix C_t = NULL ; // transpose of C
    GrB_Vector x_count = NULL;
    GxB_Iterator it = NULL ;

    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    uint64_t queue_head = 0, queue_tail = n, queue_size = n + 1;

    // TODO: decide if queue should be made inside this function.
    GRB_TRY (GrB_Matrix_new (&C_t, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_transpose (C_t, NULL, NULL, C, NULL)) ;

    // Find number of in-cluster connections for each node in the cluster
    GRB_TRY (GrB_mxm (C_t, C_t, GrB_PLUS_FP64, GxB_PLUS_SECOND_FP64, C_t, A, GrB_DESC_ST1)) ;

    // Identify nodes with in-cluster connections (more than 0)
    // TODO: use a smarter selection using c_deg
    GRB_TRY (GrB_Vector_new (&x_count, GrB_FP64, n));
    GRB_TRY (GrB_assign (x_count, NULL, NULL, 0.0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_reduce (x_count, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, C_t, GrB_DESC_T0)) ;

    uint64_t n_communities;

    GRB_TRY (GrB_Vector_new (&k, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&s_deg, GrB_FP64, n)) ;
    GRB_TRY (GrB_Matrix_new (&X, GrB_FP64, n, n)) ;
    GRB_TRY (GxB_Iterator_new(&it)) ;

    GRB_TRY (GrB_select (C_t, NULL, NULL, GrB_VALUEGT_FP64, C_t, 0.0, NULL)) ;

    // NOTE: I am wiring this so that it will be relatively easy to parallelize
    // this phase in the FUTURE

    while (queue_head != queue_tail) { // while queue not empty
        double count;
        node_id = queue[queue_head];
        queue_head = (queue_head + 1) % queue_size;
        GRB_TRY (GrB_Vector_extractElement_FP64 (&count, x_count, node_id));
        if (count <= 0.0) continue; // skip singletons

        double com_deg, n_deg;
        com_id = community [node_id] ;
        uint64_t sub_com_id = sub_com [node_id] ;

        // TODO: this is probably much easier with extract and if we select
        // nodes via their communities.
        GRB_TRY (GrB_Matrix_clear (X));
        GRB_TRY (GrB_Matrix_setElement_FP64 (X, 0.0, com_id, node_id));
        GRB_TRY (GrB_mxm (X, C_t, GrB_PLUS_FP64, GxB_PLUS_SECOND_FP64, X, A, GrB_DESC_S)) ;
        // give X the number of edges connecting x to community i
        GRB_TRY (GrB_mxm (X, NULL, NULL, GxB_PLUS_SECOND_FP64, X, S, NULL)) ;

        // momentarily remove node degree from the community total
        GRB_TRY (GrB_Vector_extractElement_FP64 (&n_deg, deg, node_id));

        // s_deg [sub_com_id] -= ndeg
        GRB_TRY (GrB_assign (s_deg, NULL, GrB_MINUS_FP64, n_deg, &sub_com_id, 1, NULL)) ;

        GRB_TRY (GrB_Vector_setElement_FP64(k, m_inv2 * n_deg, sub_com_id)) ;
        // Calculate gain in each neighboring community
        GRB_TRY (GrB_mxm (X, X, GrB_PLUS_FP64, GrB_PLUS_TIMES_SEMIRING_FP64,
            (GrB_Matrix) s_deg, (GrB_Matrix) k, GrB_DESC_ST1)) ;

        uint64_t new_c = com_id, row = 0;

        // FUTURE: This part will have to change if we want parrallel clusters
        // to work

        // Select from X (those with gain >= current and choose randomly)
        double current_gain, total_gain;
        GRB_TRY (GrB_Matrix_extractElement_FP64 (&current_gain, X, com_id, node_id)) ;
        GRB_TRY (GrB_apply (X, NULL, NULL, GrB_MINUS_FP64, X, current_gain, NULL)) ;
        GRB_TRY (GrB_select (X, NULL, NULL, GrB_VALUEGE_FP64, X, 0.0, NULL)) ;
        GRB_TRY (GrB_apply (X, NULL, NULL, GxB_POW_FP64, M_E, X, NULL)) ;
        GRB_TRY (GrB_reduce (&total_gain, NULL, GrB_PLUS_MONOID_FP64, X, NULL)) ;
        double rand_f = drand48 () * total_gain ;

        GRB_TRY (GxB_Matrix_Iterator_attach(it, X, NULL)) ;
        GrB_Info info = GxB_Matrix_Iterator_seek (it, 0);
        total_gain = 0;
        while (info == GrB_SUCCESS) {
            total_gain += GxB_Iterator_get_FP64 (it) ;

            if (total_gain > rand_f) {
                GxB_Matrix_Iterator_getIndex(it, &row, &new_c) ;
                break;
            }

            info = GxB_Matrix_Iterator_next (it);
        }

        sub_com [node_id] = new_c;
        // s_deg [new_c] += ndeg
        GRB_TRY (GrB_assign (s_deg, NULL, GrB_PLUS_FP64, n_deg, &new_c, 1, NULL));
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
    GrB_free (&C);                                          \
    GrB_free (&S);                                          \
    if (free_A) GrB_free (&A) ;                             \
    GrB_free (&mask);                                       \
    GrB_free (&deg);                                        \
    GrB_free (&c_deg);                                      \
    GrB_free (&gain);                                       \
    GrB_free (&x);                                          \
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
    LAGraph_Graph G,        // input graph (must be symmetric, no self-loops)
    uint64_t seed,          // random seed (reserved; not yet used)
    char *msg
)
{
#ifdef LG_SUITESPARSE_GRAPHBLAS_V10
    uint64_t *queue = NULL; // circular buffer with next nodes to look up
    bool *enqueued = NULL; // true if node is in the queue
    uint64_t *community = NULL ; // array: i in community c[i]
    GrB_Matrix C = NULL ; // contains readonly pointer to community array
    GrB_Matrix mask = NULL ; // mask each community
    GxB_Container cont = NULL;

    uint64_t *sub_com = NULL ; // array: i in sub-community sub_com[i]
    GrB_Matrix S = NULL ; // contains readonly pointer to sub-community array

    GrB_Vector deg = NULL ; // degree of each node
    GrB_Vector c_deg = NULL ; // degree of each community

    GrB_Matrix A = NULL ; // G->A
    GrB_Vector gain = NULL ; // gain from moving into a community
    GrB_Vector x = NULL; // full boolean vector
    GrB_Matrix node_map = NULL ; // maps aggregate node to nodes it contains
    GrB_Matrix new_node_map = NULL ; // maps aggregate node to nodes it contains

    uint64_t n, original_n;

    GxB_Iterator neighbor_it = NULL;

    bool free_A = false ;
    A = G->A;
    GRB_TRY (GrB_Matrix_nrows (&n, A));
    original_n = n;  // Store original number of nodes

    // Input checking
    LG_ASSERT_MSG (G->out_degree != NULL, LAGRAPH_NOT_CACHED,
                   "G->out_degree must be defined") ;
    LG_ASSERT_MSG(
        G->is_symmetric_structure == LAGraph_TRUE,
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be symmetric") ;

    // initialize queue
    LG_TRY (LAGraph_Malloc ((void **) &enqueued, n, sizeof(bool), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &queue, n + 1, sizeof (uint64_t), msg)) ;

    // Initialize GraphBLAS data structures
    GRB_TRY (GxB_Iterator_new(&neighbor_it));
    GRB_TRY (GxB_Container_new (&cont)) ;
    GRB_TRY (GrB_Vector_new(&x, GrB_FP64, n));
    GRB_TRY (GrB_Vector_new(&deg, GrB_FP64, n));
    GRB_TRY (GrB_Vector_new(&gain, GrB_FP64, n));

    GRB_TRY (GrB_Vector_new(&x, GrB_BOOL, n));
    GRB_TRY (GrB_assign (x, NULL, NULL, (bool) 0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_Matrix_diag (&C, x, 0)) ;
    GRB_TRY (GrB_Matrix_dup (&S, C)) ;

    // Initialize full degree vector from graph
    GRB_TRY (GrB_assign (deg, NULL, NULL, 0.0, GrB_ALL, n, NULL));
    GRB_TRY (GrB_assign (deg, NULL, GrB_PLUS_FP64, G->out_degree, GrB_ALL, n, NULL));


    double m = 0;
    GRB_TRY (GrB_reduce (&m, NULL, GrB_PLUS_MONOID_FP64, deg, NULL)) ;
    double m_inv2 = -1 / (2 * m);

    srand(seed);

    // Outer loop: repeat phases until convergence
    double prev_modularity = 0;
    double modularity_epsilon = 1e-6;

    for (int count = 0; count < LEIDEN_MAX_ITER; count++) {
        printf ("Iteration %d\n", count) ;
        // Initialize queue in random order (seed-based)
        uint64_t queue_len = n;
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

        // get i vector from C and S to use as community and sub_com arrays
        // keep inside C and S as read-only
        do {
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
            GRB_TRY (GxB_unload_Matrix_into_Container (S, cont, NULL)) ;
            GRB_TRY (GxB_Vector_unload (cont->i, (void **) &sub_com, &type, &n_community,
                &X_memsize, &handling, NULL)) ;
            GRB_TRY (GxB_Vector_load (cont->i, (void **) &sub_com, type, n_community,
                X_memsize, GxB_IS_READONLY, NULL)) ;
            GRB_TRY (GxB_load_Matrix_from_Container(S, cont, NULL)) ;
        } while (0) ;

        //----------------------------------------------------------------------
        // Phase 1: move nodes
        //----------------------------------------------------------------------
        LG_TRY (LG_Leiden_move_nodes (C, community, A, deg, queue, enqueued, m_inv2, msg)) ;

        //----------------------------------------------------------------------
        // Prep phase 2
        //----------------------------------------------------------------------
        // Initialize queue
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

        //----------------------------------------------------------------------
        // Phase 2: refine
        //----------------------------------------------------------------------
        // Save Phase 1 communities as parent communities
        LG_TRY (LG_Leiden_refinement (S, sub_com, C, community, c_deg, A, deg, queue, m_inv2, msg)) ;

        //----------------------------------------------------------------------
        // Phase 3: aggregate
        //----------------------------------------------------------------------

        // Compute modularity to check convergence
        double current_modularity;
        LG_TRY (LAGr_AdjModularity (&current_modularity, 1.0, A, S, msg));

        ASSERT (current_modularity >= prev_modularity) ;
        if (fabs(current_modularity - prev_modularity) < modularity_epsilon) {
            break;  // converged
        }
        prev_modularity = current_modularity;

        LG_TRY (LG_Leiden_aggregate (&S, &C, &A, free_A, msg)) ;
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
        LG_TRY (LAGr_AdjModularity (&current_modularity, 1.0, A, S, msg));
        ASSERT (current_modularity == prev_modularity) ;
    #endif
    } // end outer while loop


    // Map refined communities back to original nodes
    GRB_TRY (GrB_Vector_new(c_handle, GrB_INT64, original_n));
    GRB_TRY (GrB_free (&C)) ;
    if (node_map == NULL) {
        ASSERT (n == original_n) ;
        GRB_TRY (GxB_Vector_load (*c_handle, (void **) &community, GrB_INT64, n,
            sizeof (int64_t) * n, GrB_DEFAULT, NULL));
    } else {
        GRB_TRY (GxB_Matrix_Iterator_attach (neighbor_it, node_map, NULL)) ;
        GrB_Info info = GxB_Matrix_Iterator_seek (neighbor_it, 0) ;
        while (info == GrB_SUCCESS) {
            uint64_t coarse_id, orig_id;
            GxB_Matrix_Iterator_getIndex (neighbor_it, &coarse_id, &orig_id) ;
            GRB_TRY (GrB_Vector_setElement_INT64 (*c_handle, community[coarse_id], orig_id)) ;
            info = GxB_Matrix_Iterator_next (neighbor_it) ;
        }
    }

    LG_FREE_WORK;
    return GrB_SUCCESS;
#else
    return GrB_NOT_IMPLEMENTED;
#endif
}

