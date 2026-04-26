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
//       score(i->c) = T[c] - k[i] * k_comm[c] / m
//     where T[c] = sum_{j in c, j!=i} A[i,j], k_comm[c] = total degree of
//     community c (excluding i), and m = total edge weight / 2 (constant).
//     Initial partition is induced by the Phase-1 communities of the previous
//     level (singletons on the first level).
//     Repeat until no node changes community.
//
//   Phase 2 (Refinement – key Leiden addition): Save Phase-1 communities as
//     "parent" communities.  Restart each node in its own singleton
//     sub-community.  In each local-move step a node may only join a
//     sub-community whose parent equals its own Phase-1 parent.  This ensures
//     every output community is an internally well-connected subgraph of the
//     corresponding Phase-1 community.
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
    GrB_free (&k_vec) ;                                     \
    GrB_free (&v) ;                                         \
    GrB_free (&A_agg) ;                                     \
    GrB_free (&S_mat) ;                                     \
    GrB_free (&A_temp) ;                                    \
    LAGraph_Free ((void **) &k_arr,      NULL) ;            \
    LAGraph_Free ((void **) &c_arr,      NULL) ;            \
    LAGraph_Free ((void **) &k_comm,     NULL) ;            \
    LAGraph_Free ((void **) &c_p1,       NULL) ;            \
    LAGraph_Free ((void **) &c_ref,      NULL) ;            \
    LAGraph_Free ((void **) &k_ref_comm, NULL) ;            \
    LAGraph_Free ((void **) &T_local,    NULL) ;            \
    LAGraph_Free ((void **) &dirty,      NULL) ;            \
    LAGraph_Free ((void **) &dirty_list, NULL) ;            \
    LAGraph_Free ((void **) &nbrs_j,     NULL) ;            \
    LAGraph_Free ((void **) &nbrs_v,     NULL) ;            \
    LAGraph_Free ((void **) &remap,      NULL) ;            \
    LAGraph_Free ((void **) &o_comm,     NULL) ;            \
    LAGraph_Free ((void **) &init_comm,  NULL) ;            \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
    if (c_handle != NULL) GrB_free (c_handle) ;             \
}

#include "LG_internal.h"
#include <LAGraphX.h>
#include <LAGraph.h>
#include <stdlib.h>
#include <string.h>

#define LEIDEN_MAX_ITER 100

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

    //--------------------------------------------------------------------------
    // declare all workspace (must precede any LG_TRY/GRB_TRY calls so that
    // LG_FREE_ALL can safely free them even on early exit)
    //--------------------------------------------------------------------------

    GrB_Vector  k_vec      = NULL ;
    GrB_Vector  v          = NULL ;
    GrB_Matrix  A_agg      = NULL ;   // owned coarsened graph (Phase 3)
    GrB_Matrix  S_mat      = NULL ;   // temporary membership matrix (Phase 3)
    GrB_Matrix  A_temp     = NULL ;   // temporary for mxm (Phase 3)
    double     *k_arr      = NULL ;   // k_arr[i]      = degree of node i (current level)
    int64_t    *c_arr      = NULL ;   // c_arr[i]      = Phase-1 community label
    double     *k_comm     = NULL ;   // k_comm[l]     = total degree of community l
    int64_t    *c_p1       = NULL ;   // c_p1[i]       = Phase-1 parent community
    int64_t    *c_ref      = NULL ;   // c_ref[i]      = refined sub-community label
    double     *k_ref_comm = NULL ;   // k_ref_comm[l] = total degree of sub-community l
    double     *T_local    = NULL ;   // scratch: edge sums from node i to each community
    int8_t     *dirty      = NULL ;   // dirty[l] = 1 if T_local[l] was written
    GrB_Index  *dirty_list = NULL ;   // list of community labels touched this node
    GrB_Index  *nbrs_j     = NULL ;   // scratch: extracted neighbor indices
    double     *nbrs_v     = NULL ;   // scratch: extracted neighbor weights
    GrB_Index  *remap      = NULL ;   // remap[old_label] -> new contiguous label
    GrB_Index  *o_comm     = NULL ;   // o_comm[i] = community of original node i
    GrB_Index  *init_comm  = NULL ;   // init_comm[r] = initial c_arr for aggregate node r

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (c_handle != NULL, GrB_NULL_POINTER) ;
    (*c_handle) = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        (G->kind == LAGraph_ADJACENCY_DIRECTED &&
         G->is_symmetric_structure == LAGraph_TRUE),
        LAGRAPH_NOT_CACHED,
        "G must be undirected or have symmetric structure") ;

    GrB_Matrix A = G->A ;
    GrB_Index  n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    // Degenerate: return immediately for empty (0-node) graph.
    if (n == 0)
    {
        GRB_TRY (GrB_Vector_new (c_handle, GrB_INT64, 0)) ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // allocate workspace (all arrays sized n; used for indices 0..n_cur-1)
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Malloc ((void **) &k_arr,      n, sizeof (double),    msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &c_arr,      n, sizeof (int64_t),   msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &k_comm,     n, sizeof (double),    msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &c_p1,       n, sizeof (int64_t),   msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &c_ref,      n, sizeof (int64_t),   msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &k_ref_comm, n, sizeof (double),    msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &T_local,    n, sizeof (double),    msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &dirty,      n, sizeof (int8_t),    msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &dirty_list, n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &nbrs_j,     n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &nbrs_v,     n, sizeof (double),    msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &remap,      n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &o_comm,     n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &init_comm,  n, sizeof (GrB_Index), msg)) ;

    //--------------------------------------------------------------------------
    // compute m = total edge weight / 2 from G->A (invariant under aggregation)
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&k_vec, GrB_FP64, n)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (k_vec, NULL, NULL,
        GrB_PLUS_MONOID_FP64, A, NULL)) ;
    double m = 0.0 ;
    GRB_TRY (GrB_Vector_reduce_FP64 (&m, NULL, GrB_PLUS_MONOID_FP64,
        k_vec, NULL)) ;
    m /= 2.0 ;
    GrB_free (&k_vec) ;
    k_vec = NULL ;

    // Empty graph: return a singleton partition.
    if (m == 0.0)
    {
        GRB_TRY (GrB_Vector_new (c_handle, GrB_INT64, n)) ;
        for (GrB_Index i = 0 ; i < n ; i++)
        {
            GRB_TRY (GrB_Vector_setElement_INT64 (*c_handle, (int64_t) i, i)) ;
        }
        LG_FREE_WORK ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // initialise multi-level state
    //
    // o_comm[i]   = community of original node i (tracks the composition of
    //               all levels' refined partitions).  Starts as identity.
    // init_comm[r] = initial Phase-1 community for aggregate node r at the
    //               next level.  First level: singleton start (c_arr[i] = i).
    //--------------------------------------------------------------------------

    for (GrB_Index i = 0 ; i < n ; i++)
    {
        o_comm[i]    = i ;
        init_comm[i] = i ;
    }

    GrB_Matrix A_cur = A ;   // current-level graph (not owned: either G->A or A_agg)
    GrB_Index  n_cur = n ;

    //==========================================================================
    // OUTER AGGREGATION LOOP
    //==========================================================================

    bool outer_changed = true ;
    while (outer_changed)
    {
        outer_changed = false ;

        //----------------------------------------------------------------------
        // create GraphBLAS vectors sized for the current level
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Vector_new (&k_vec, GrB_FP64, n_cur)) ;
        GRB_TRY (GrB_Vector_new (&v,     GrB_FP64, n_cur)) ;

        //----------------------------------------------------------------------
        // compute k_arr[0..n_cur-1] from A_cur (includes self-loops in A_agg)
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Matrix_reduce_Monoid (k_vec, NULL, NULL,
            GrB_PLUS_MONOID_FP64, A_cur, NULL)) ;
        for (GrB_Index i = 0 ; i < n_cur ; i++)
        {
            GrB_Info info = GrB_Vector_extractElement_FP64 (&k_arr[i], k_vec, i) ;
            if (info == GrB_NO_VALUE) k_arr[i] = 0.0 ;
        }

        //----------------------------------------------------------------------
        // PHASE 1: Local Move Phase
        //
        // Initialise partition from init_comm (singletons on first level;
        // induced Phase-1 partition on subsequent levels).
        // Score: score(i->c) = T[c] - k[i]*k_comm[c]/m  (self-loops skipped).
        //----------------------------------------------------------------------

        memset (dirty,   0, n * sizeof (int8_t)) ;
        memset (T_local, 0, n * sizeof (double)) ;
        memset (k_comm,  0, n * sizeof (double)) ;
        for (GrB_Index i = 0 ; i < n_cur ; i++)
        {
            c_arr[i]          = (int64_t) init_comm[i] ;
            k_comm[init_comm[i]] += k_arr[i] ;
        }

        bool changed = true ;
        for (int p1_iter = 0 ; changed && p1_iter < LEIDEN_MAX_ITER ; p1_iter++)
        {
            changed = false ;
            for (GrB_Index i = 0 ; i < n_cur ; i++)
            {
                double ki = k_arr[i] ;
                if (ki == 0.0) continue ;

                int64_t ci = c_arr[i] ;

                GRB_TRY (GrB_Col_extract (v, NULL, NULL, A_cur,
                    GrB_ALL, n_cur, i, GrB_DESC_T0)) ;

                GrB_Index nvals ;
                GRB_TRY (GrB_Vector_nvals (&nvals, v)) ;
                if (nvals == 0) continue ;

                GRB_TRY (GrB_Vector_extractTuples_FP64 (
                    nbrs_j, nbrs_v, &nvals, v)) ;

                // Temporarily remove i from community ci.
                k_comm[ci] -= ki ;

                GrB_Index ndirty = 0 ;
                for (GrB_Index t = 0 ; t < nvals ; t++)
                {
                    GrB_Index j = nbrs_j[t] ;
                    if (j == i) continue ;          // skip self-loop (in A_agg)
                    int64_t cj = c_arr[j] ;
                    if (!dirty[cj])
                    {
                        dirty[cj]            = 1 ;
                        dirty_list[ndirty++] = (GrB_Index) cj ;
                        T_local[cj]          = 0.0 ;
                    }
                    T_local[cj] += nbrs_v[t] ;
                }

                double  T_ci      = dirty[ci] ? T_local[ci] : 0.0 ;
                double  score_ci  = T_ci - ki * k_comm[ci] / m ;
                double  best_score = score_ci ;
                int64_t best_c     = ci ;

                for (GrB_Index d = 0 ; d < ndirty ; d++)
                {
                    int64_t c_cand = (int64_t) dirty_list[d] ;
                    if (c_cand == ci) continue ;
                    double score = T_local[c_cand] - ki * k_comm[c_cand] / m ;
                    if (score > best_score)
                    {
                        best_score = score ;
                        best_c     = c_cand ;
                    }
                }

                c_arr[i] = best_c ;
                if (best_c == ci)
                {
                    k_comm[ci] += ki ;
                }
                else
                {
                    k_comm[best_c] += ki ;
                    changed = true ;
                }

                for (GrB_Index d = 0 ; d < ndirty ; d++)
                {
                    dirty[dirty_list[d]] = 0 ;
                }
            }
        }

        //----------------------------------------------------------------------
        // PHASE 2: Refinement Phase (key Leiden addition)
        //
        // Save Phase-1 result.  Restart each node in a singleton sub-community.
        // Only allow moves within the same Phase-1 parent community.
        //----------------------------------------------------------------------

        memcpy (c_p1, c_arr, n_cur * sizeof (int64_t)) ;

        for (GrB_Index i = 0 ; i < n_cur ; i++)
        {
            c_ref[i]      = (int64_t) i ;
            k_ref_comm[i] = k_arr[i] ;
        }

        changed = true ;
        for (int p2_iter = 0 ; changed && p2_iter < LEIDEN_MAX_ITER ; p2_iter++)
        {
            changed = false ;
            for (GrB_Index i = 0 ; i < n_cur ; i++)
            {
                double ki = k_arr[i] ;
                if (ki == 0.0) continue ;

                int64_t pi     = c_p1[i] ;
                int64_t ci_ref = c_ref[i] ;

                GRB_TRY (GrB_Col_extract (v, NULL, NULL, A_cur,
                    GrB_ALL, n_cur, i, GrB_DESC_T0)) ;

                GrB_Index nvals ;
                GRB_TRY (GrB_Vector_nvals (&nvals, v)) ;
                if (nvals == 0) continue ;

                GRB_TRY (GrB_Vector_extractTuples_FP64 (
                    nbrs_j, nbrs_v, &nvals, v)) ;

                k_ref_comm[ci_ref] -= ki ;

                GrB_Index ndirty = 0 ;
                for (GrB_Index t = 0 ; t < nvals ; t++)
                {
                    GrB_Index j = nbrs_j[t] ;
                    if (j == i) continue ;              // skip self-loop
                    if (c_p1[j] != pi) continue ;       // cross-parent: skip

                    int64_t cj_ref = c_ref[j] ;
                    if (!dirty[cj_ref])
                    {
                        dirty[cj_ref]        = 1 ;
                        dirty_list[ndirty++] = (GrB_Index) cj_ref ;
                        T_local[cj_ref]      = 0.0 ;
                    }
                    T_local[cj_ref] += nbrs_v[t] ;
                }

                double  T_ci_ref    = dirty[ci_ref] ? T_local[ci_ref] : 0.0 ;
                double  score_ci_ref = T_ci_ref - ki * k_ref_comm[ci_ref] / m ;
                double  best_score  = score_ci_ref ;
                int64_t best_c_ref  = ci_ref ;

                for (GrB_Index d = 0 ; d < ndirty ; d++)
                {
                    int64_t c_cand = (int64_t) dirty_list[d] ;
                    if (c_cand == ci_ref) continue ;
                    double score = T_local[c_cand] - ki * k_ref_comm[c_cand] / m ;
                    if (score > best_score)
                    {
                        best_score = score ;
                        best_c_ref = c_cand ;
                    }
                }

                c_ref[i] = best_c_ref ;
                if (best_c_ref == ci_ref)
                {
                    k_ref_comm[ci_ref] += ki ;
                }
                else
                {
                    k_ref_comm[best_c_ref] += ki ;
                    changed = true ;
                }

                for (GrB_Index d = 0 ; d < ndirty ; d++)
                {
                    dirty[dirty_list[d]] = 0 ;
                }
            }
        }

        //----------------------------------------------------------------------
        // Relabel c_ref to contiguous integers 0..K_ref-1
        //----------------------------------------------------------------------

        // Use n as sentinel ("not yet assigned"); safe because c_ref values
        // are in 0..n_cur-1 < n.
        for (GrB_Index i = 0 ; i < n_cur ; i++) remap[i] = n ;

        GrB_Index K_ref = 0 ;
        for (GrB_Index i = 0 ; i < n_cur ; i++)
        {
            GrB_Index old_label = (GrB_Index) c_ref[i] ;
            if (remap[old_label] == n) remap[old_label] = K_ref++ ;
            c_ref[i] = (int64_t) remap[old_label] ;
        }

        //----------------------------------------------------------------------
        // Compute init_comm for next level.
        //
        // Aggregate node r (0..K_ref-1) is the refined community c_ref[i] for
        // any i with that label.  All such nodes have the same Phase-1 parent
        // c_arr[i] (Leiden invariant), so we record that as the initial
        // community for aggregate node r in the next outer iteration.
        //----------------------------------------------------------------------

        for (GrB_Index i = 0 ; i < n_cur ; i++)
        {
            // c_ref[i] is in 0..K_ref-1 and init_comm is size n >= K_ref.
            init_comm[(GrB_Index) c_ref[i]] = (GrB_Index) c_arr[i] ;
        }

        //----------------------------------------------------------------------
        // Compose o_comm: original node i now maps to aggregate community
        // c_ref[o_comm[i]].  o_comm[i] is always a valid index in c_ref
        // because it was set to some value in 0..n_cur-1 in the previous
        // iteration (or to i on the first iteration).
        //----------------------------------------------------------------------

        for (GrB_Index i = 0 ; i < n ; i++)
        {
            o_comm[i] = (GrB_Index) c_ref[o_comm[i]] ;
        }

        //----------------------------------------------------------------------
        // PHASE 3: Aggregation — build coarsened graph if communities merged
        //----------------------------------------------------------------------

        GrB_free (&k_vec) ;  k_vec = NULL ;
        GrB_free (&v) ;      v     = NULL ;

        if (K_ref < n_cur)
        {
            outer_changed = true ;

            // S_mat: n_cur × K_ref indicator matrix; S[i, c_ref[i]] = 1
            GRB_TRY (GrB_Matrix_new (&S_mat, GrB_FP64, n_cur, K_ref)) ;
            for (GrB_Index i = 0 ; i < n_cur ; i++)
            {
                GRB_TRY (GrB_Matrix_setElement_FP64 (
                    S_mat, 1.0, i, (GrB_Index) c_ref[i])) ;
            }
            GRB_TRY (GrB_Matrix_wait (S_mat, GrB_MATERIALIZE)) ;

            // A_temp = A_cur * S  (n_cur × K_ref)
            GRB_TRY (GrB_Matrix_new (&A_temp, GrB_FP64, n_cur, K_ref)) ;
            GRB_TRY (GrB_mxm (A_temp, NULL, NULL,
                GrB_PLUS_TIMES_SEMIRING_FP64, A_cur, S_mat, NULL)) ;

            // A_new = S^T * A_temp  (K_ref × K_ref)
            GrB_Matrix A_new = NULL ;
            GRB_TRY (GrB_Matrix_new (&A_new, GrB_FP64, K_ref, K_ref)) ;
            GRB_TRY (GrB_mxm (A_new, NULL, NULL,
                GrB_PLUS_TIMES_SEMIRING_FP64, S_mat, A_temp, GrB_DESC_T0)) ;

            GrB_free (&S_mat) ;  S_mat  = NULL ;
            GrB_free (&A_temp) ; A_temp = NULL ;
            GrB_free (&A_agg) ;  // free previous level's aggregate graph
            A_agg  = A_new ;
            A_cur  = A_agg ;
            n_cur  = K_ref ;
        }
        // K_ref == n_cur: no communities merged this level → converged.
    }

    //--------------------------------------------------------------------------
    // Build output GrB_Vector from o_comm
    // (o_comm values are already relabeled 0..K_final-1 from the last iteration)
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (c_handle, GrB_INT64, n)) ;
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        GRB_TRY (GrB_Vector_setElement_INT64 (*c_handle, (int64_t) o_comm[i], i)) ;
    }

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
