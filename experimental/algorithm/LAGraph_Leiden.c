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

#include "LG_internal.h"
#include <LAGraphX.h>
#include <LAGraph.h>
#include <stdlib.h>
#include <string.h>

#define LEIDEN_MAX_ITER 100

// The CSR fast path uses the SuiteSparse:GraphBLAS Container API
// (GxB_Container, GxB_unload_Matrix_into_Container, GxB_Vector_load /
// GxB_Vector_unload, GxB_load_Matrix_from_Container), introduced in
// SuiteSparse:GraphBLAS v10.0.0.  On older versions we fall back to a
// CSR materialization via GrB_Matrix_extractTuples + counting-sort scatter.
#ifndef LAGR_LEIDEN_USE_CONTAINER
#if defined(GxB_IMPLEMENTATION) && (GxB_IMPLEMENTATION >= GxB_VERSION(10,0,0))
#define LAGR_LEIDEN_USE_CONTAINER 1
#else
#define LAGR_LEIDEN_USE_CONTAINER 0
#endif
#endif

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
    GrB_free (&k_vec) ;                                     \
    GrB_free (&A_agg) ;                                     \
    GrB_free (&A_new) ;                                     \
    GrB_free (&S_mat) ;                                     \
    GrB_free (&A_temp) ;                                    \
    GrB_free (&one_scalar) ;                                \
    LAGR_LEIDEN_FREE_CONTAINER ;                            \
    LAGraph_Free ((void **) &k_arr,      NULL) ;            \
    LAGraph_Free ((void **) &c_arr,      NULL) ;            \
    LAGraph_Free ((void **) &k_comm,     NULL) ;            \
    LAGraph_Free ((void **) &c_p1,       NULL) ;            \
    LAGraph_Free ((void **) &c_ref,      NULL) ;            \
    LAGraph_Free ((void **) &k_ref_comm, NULL) ;            \
    LAGraph_Free ((void **) &T_local,    NULL) ;            \
    LAGraph_Free ((void **) &dirty,      NULL) ;            \
    LAGraph_Free ((void **) &dirty_list, NULL) ;            \
    LAGraph_Free ((void **) &remap,      NULL) ;            \
    LAGraph_Free ((void **) &o_comm,     NULL) ;            \
    LAGraph_Free ((void **) &init_comm,  NULL) ;            \
    LAGraph_Free ((void **) &Ap,         NULL) ;            \
    LAGraph_Free ((void **) &Aj,         NULL) ;            \
    LAGraph_Free ((void **) &Ax,         NULL) ;            \
    LAGraph_Free ((void **) &I_tup,      NULL) ;            \
    LAGraph_Free ((void **) &J_tup,      NULL) ;            \
    LAGraph_Free ((void **) &X_tup,      NULL) ;            \
    LAGraph_Free ((void **) &cursor,     NULL) ;            \
    LAGraph_Free ((void **) &iota,       NULL) ;            \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
    if (c_handle != NULL) GrB_free (c_handle) ;             \
}

#if LAGR_LEIDEN_USE_CONTAINER
#define LAGR_LEIDEN_FREE_CONTAINER GxB_Container_free (&cont)
#else
#define LAGR_LEIDEN_FREE_CONTAINER ((void) 0)
#endif

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

    GrB_Vector    k_vec      = NULL ;
    GrB_Matrix    A_agg      = NULL ;   // owned coarsened graph (Phase 3)
    GrB_Matrix    A_new      = NULL ;   // next-level aggregate before ownership transfer
    GrB_Matrix    S_mat      = NULL ;   // temporary membership matrix (Phase 3)
    GrB_Matrix    A_temp     = NULL ;   // temporary for mxm (Phase 3)
    GrB_Scalar    one_scalar = NULL ;   // FP64 scalar with value 1.0 for build_Scalar
#if LAGR_LEIDEN_USE_CONTAINER
    GxB_Container cont       = NULL ;   // for unloading A_cur into raw CSR arrays
#endif
    double       *k_arr      = NULL ;   // k_arr[i]      = degree of node i (current level)
    int64_t      *c_arr      = NULL ;   // c_arr[i]      = Phase-1 community label
    double       *k_comm     = NULL ;   // k_comm[l]     = total degree of community l
    int64_t      *c_p1       = NULL ;   // c_p1[i]       = Phase-1 parent community
    int64_t      *c_ref      = NULL ;   // c_ref[i]      = refined sub-community label
    double       *k_ref_comm = NULL ;   // k_ref_comm[l] = total degree of sub-community l
    double       *T_local    = NULL ;   // scratch: edge sums from node i to each community
    int8_t       *dirty      = NULL ;   // dirty[l] = 1 if T_local[l] was written
    GrB_Index    *dirty_list = NULL ;   // list of community labels touched this node
    GrB_Index    *remap      = NULL ;   // remap[old_label] -> new contiguous label
    int64_t      *o_comm     = NULL ;   // o_comm[i] = community of original node i
    GrB_Index    *init_comm  = NULL ;   // init_comm[r] = initial c_arr for aggregate node r

    // Raw CSR pointers for inner-loop walks.  On v10+ they are obtained by
    // unloading A_cur into the SuiteSparse Container (zero-copy); ownership
    // returns to GraphBLAS on reload (we then null them so LG_FREE_WORK
    // doesn't double-free).  On older versions they are allocated by us
    // and (re)filled per level via GrB_Matrix_extractTuples + counting-sort.
    GrB_Index  *Ap         = NULL ;   // row pointers, size n_cur+1
    GrB_Index  *Aj         = NULL ;   // column indices, size Anz
    double     *Ax         = NULL ;   // values, size Anz (only if !iso on v10)
    GrB_Index  *I_tup      = NULL ;   // raw row indices from extractTuples (v9 fallback)
    GrB_Index  *J_tup      = NULL ;   // raw col indices from extractTuples (v9 fallback)
    double     *X_tup      = NULL ;   // raw values from extractTuples (v9 fallback)
    GrB_Index  *cursor     = NULL ;   // scatter cursor for CSR build (v9 fallback)
    GrB_Index  *iota       = NULL ;   // [0,1,...,n-1] for vector/matrix build
    GrB_Index   Ap_cap     = 0 ;      // current allocated capacity of Ap (v9 fallback)
    GrB_Index   Anz_cap    = 0 ;      // current allocated capacity of Aj/Ax/tuples (v9 fallback)

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    (void) seed ;       // reserved for future randomized refinement
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
    LG_TRY (LAGraph_Malloc ((void **) &remap,      n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &o_comm,     n, sizeof (int64_t),   msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &init_comm,  n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &iota,       n, sizeof (GrB_Index), msg)) ;
    for (GrB_Index i = 0 ; i < n ; i++) iota[i] = i ;

    // Reusable FP64 scalar with value 1.0 for GxB_Matrix_build_Scalar.
    GRB_TRY (GrB_Scalar_new (&one_scalar, GrB_FP64)) ;
    GRB_TRY (GrB_Scalar_setElement_FP64 (one_scalar, 1.0)) ;

    // Reusable container for unloading A_cur into raw CSR arrays per level
    // (only available with the SuiteSparse v10+ Container API).
#if LAGR_LEIDEN_USE_CONTAINER
    GRB_TRY (GxB_Container_new (&cont)) ;
#endif

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
    double two_m = 2.0 * m ;        // denominator of the modularity penalty
    GrB_free (&k_vec) ;
    k_vec = NULL ;

    // Empty graph: return a singleton partition.
    if (m == 0.0)
    {
        // c[i] = i for all i, built in one call instead of n setElement calls.
        for (GrB_Index i = 0 ; i < n ; i++) c_arr[i] = (int64_t) i ;
        GRB_TRY (GrB_Vector_new (c_handle, GrB_INT64, n)) ;
        GRB_TRY (GrB_Vector_build_INT64 (*c_handle, iota, c_arr, n,
            GrB_FIRST_INT64)) ;
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
        o_comm[i]    = (int64_t) i ;
        init_comm[i] = i ;
    }

    // Duplicate G->A as FP64 so we own A_cur and may unload it via the
    // container API expecting FP64 values.  Done after the m == 0 early
    // return to avoid an unused copy for empty-edge graphs.  G->A may be
    // any numeric type (BOOL on pattern-only matrices, INT*, FP32, ...);
    // we typecast once here so the inner-loop CSR walks always read double.
    GRB_TRY (GrB_Matrix_new (&A_agg, GrB_FP64, n, n)) ;
    GRB_TRY (GrB_Matrix_assign (A_agg, NULL, NULL, A,
        GrB_ALL, n, GrB_ALL, n, NULL)) ;
    GrB_Matrix A_cur = A_agg ;
    GrB_Index  n_cur = n ;

    //==========================================================================
    // OUTER AGGREGATION LOOP
    //==========================================================================

    bool outer_changed = true ;
    while (outer_changed)
    {
        outer_changed = false ;

        //----------------------------------------------------------------------
        // Compute degrees k_arr[i] = sum of row i (includes self-loops in
        // A_agg) and unload A_cur into raw CSR arrays Ap/Aj/Ax via the
        // SuiteSparse Container API.  This replaces O(n_cur) GrB_Col_extract
        // + extractTuples calls per inner iteration with direct pointer walks.
        //
        // Force A_cur into a deterministic state before unloading: dense
        // reduce target for k, sparse + row-major + non-iso + 64-bit indices
        // for A.  Container fields are then asserted to match expectations.
        //----------------------------------------------------------------------

        // 1) Compute degrees with GraphBLAS reduce *before* unloading the
        //    matrix.  Zero-fill k_vec first so isolated rows produce 0.0
        //    (otherwise reduce leaves them as missing entries).
        GrB_free (&k_vec) ;
        GRB_TRY (GrB_Vector_new (&k_vec, GrB_FP64, n_cur)) ;
        GRB_TRY (GrB_assign (k_vec, NULL, NULL, (double) 0.0,
            GrB_ALL, n_cur, NULL)) ;
        GRB_TRY (GrB_Matrix_reduce_Monoid (k_vec, NULL, GrB_PLUS_FP64,
            GrB_PLUS_MONOID_FP64, A_cur, NULL)) ;

#if LAGR_LEIDEN_USE_CONTAINER
        // Unload k_vec into k_arr (dense FP64 array of length n_cur).
        // The previous k_arr workspace allocation is replaced by a pointer
        // owned by GraphBLAS until LAGraph_Free reclaims it.
        LAGraph_Free ((void **) &k_arr, NULL) ;
        {
            GrB_Type k_type = NULL ;
            uint64_t k_n = 0, k_size = 0 ;
            int      k_handling = GrB_DEFAULT ;
            void    *k_void = NULL ;
            GRB_TRY (GxB_Vector_unload (k_vec, &k_void, &k_type, &k_n,
                &k_size, &k_handling, NULL)) ;
            LG_ASSERT_MSG (k_type == GrB_FP64 && k_n == n_cur,
                GrB_INVALID_VALUE,
                "k_vec unload: unexpected type or length") ;
            k_arr = (double *) k_void ;
        }

        // 2) Force A_cur into the format we want, then unload into container.
        //    Hints: sparse, row-major, non-iso, 64-bit row pointers/indices.
        //    GxB_unload_Matrix_into_Container materializes pending work.
        GRB_TRY (GrB_set (A_cur, GxB_SPARSE, GxB_SPARSITY_CONTROL)) ;
        GRB_TRY (GrB_set (A_cur, (int32_t) GrB_ROWMAJOR,
            GrB_STORAGE_ORIENTATION_HINT)) ;
        GRB_TRY (GrB_Matrix_set_INT32 (A_cur, false, GxB_ISO)) ;
        GRB_TRY (GrB_Matrix_set_INT32 (A_cur, 64, GxB_OFFSET_INTEGER_HINT)) ;
        GRB_TRY (GrB_Matrix_set_INT32 (A_cur, 64, GxB_ROWINDEX_INTEGER_HINT)) ;
        GRB_TRY (GrB_Matrix_set_INT32 (A_cur, 64, GxB_COLINDEX_INTEGER_HINT)) ;
        GRB_TRY (GrB_wait (A_cur, GrB_MATERIALIZE)) ;

        GRB_TRY (GxB_unload_Matrix_into_Container (A_cur, cont, NULL)) ;
        LG_ASSERT_MSG (cont->format == GxB_SPARSE,
            GrB_INVALID_VALUE, "A_cur container is not sparse CSR") ;
        LG_ASSERT_MSG (cont->orientation == GrB_ROWMAJOR,
            GrB_INVALID_VALUE, "A_cur container is not row-major") ;
        LG_ASSERT_MSG (!cont->iso,
            GrB_INVALID_VALUE, "A_cur container unexpectedly iso") ;

        // Unload row pointers, column indices, and values from the container's
        // internal vectors into raw arrays for the inner-loop CSR walks.
        // We hold these arrays as our own until reload below.
        GrB_Type   pty = NULL,  ity = NULL,  xty = NULL ;
        uint64_t   pn  = 0,     in_  = 0,    xn  = 0 ;
        uint64_t   psz = 0,     isz = 0,     xsz = 0 ;
        int        ph  = GrB_DEFAULT, ih = GrB_DEFAULT, xh = GrB_DEFAULT ;
        void      *pv  = NULL,  *iv  = NULL, *xv  = NULL ;

        GRB_TRY (GxB_Vector_unload (cont->p, &pv, &pty, &pn, &psz, &ph, NULL));
        GRB_TRY (GxB_Vector_unload (cont->i, &iv, &ity, &in_, &isz, &ih, NULL));
        GRB_TRY (GxB_Vector_unload (cont->x, &xv, &xty, &xn, &xsz, &xh, NULL));
        // Offsets must be 64-bit unsigned (we forced via INTEGER_HINT).
        // Column indices may come back as either UINT64 or INT64 depending
        // on SuiteSparse's internal choice; both have identical bit width
        // and represent non-negative indices, so reinterpret cast is safe.
        // Values must be FP64 since A_agg was constructed as FP64.
        LG_ASSERT_MSG (pty == GrB_UINT64,
            GrB_INVALID_VALUE, "container offsets are not 64-bit unsigned") ;
        LG_ASSERT_MSG (ity == GrB_UINT64 || ity == GrB_INT64,
            GrB_INVALID_VALUE, "container indices are not 64-bit") ;
        LG_ASSERT_MSG (xty == GrB_FP64,
            GrB_INVALID_VALUE, "container values are not FP64") ;
        Ap = (GrB_Index *) pv ;
        Aj = (GrB_Index *) iv ;
        Ax = (double    *) xv ;
#else
        // Fallback for SuiteSparse:GraphBLAS < v10.0.0 (no Container API):
        // copy degrees out of k_vec, then materialize CSR via extractTuples
        // + counting-sort scatter.  k_vec contains entries only for non-zero
        // rows, so zero-fill k_arr first then scatter.
        for (GrB_Index i = 0 ; i < n_cur ; i++) k_arr[i] = 0.0 ;
        {
            GrB_Index nvk = n_cur ;
            GrB_Index *Ik = NULL ;
            double    *Xk = NULL ;
            LG_TRY (LAGraph_Malloc ((void **) &Ik, n_cur, sizeof (GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void **) &Xk, n_cur, sizeof (double),    msg)) ;
            GRB_TRY (GrB_Vector_extractTuples_FP64 (Ik, Xk, &nvk, k_vec)) ;
            for (GrB_Index t = 0 ; t < nvk ; t++) k_arr[Ik[t]] = Xk[t] ;
            LAGraph_Free ((void **) &Ik, NULL) ;
            LAGraph_Free ((void **) &Xk, NULL) ;
        }

        GrB_Index Anz ;
        GRB_TRY (GrB_Matrix_nvals (&Anz, A_cur)) ;
        if (Ap_cap < n_cur + 1)
        {
            LAGraph_Free ((void **) &Ap,     NULL) ;
            LAGraph_Free ((void **) &cursor, NULL) ;
            LG_TRY (LAGraph_Malloc ((void **) &Ap,     n_cur + 1,
                sizeof (GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void **) &cursor, n_cur,
                sizeof (GrB_Index), msg)) ;
            Ap_cap = n_cur + 1 ;
        }
        if (Anz_cap < Anz)
        {
            GrB_Index newcap = (Anz < 16) ? 16 : Anz ;
            LAGraph_Free ((void **) &Aj,    NULL) ;
            LAGraph_Free ((void **) &Ax,    NULL) ;
            LAGraph_Free ((void **) &I_tup, NULL) ;
            LAGraph_Free ((void **) &J_tup, NULL) ;
            LAGraph_Free ((void **) &X_tup, NULL) ;
            LG_TRY (LAGraph_Malloc ((void **) &Aj,    newcap,
                sizeof (GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void **) &Ax,    newcap,
                sizeof (double),    msg)) ;
            LG_TRY (LAGraph_Malloc ((void **) &I_tup, newcap,
                sizeof (GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void **) &J_tup, newcap,
                sizeof (GrB_Index), msg)) ;
            LG_TRY (LAGraph_Malloc ((void **) &X_tup, newcap,
                sizeof (double),    msg)) ;
            Anz_cap = newcap ;
        }

        memset (Ap, 0, (n_cur + 1) * sizeof (GrB_Index)) ;
        if (Anz > 0)
        {
            GrB_Index nout = Anz ;
            GRB_TRY (GrB_Matrix_extractTuples_FP64 (I_tup, J_tup, X_tup,
                &nout, A_cur)) ;
            for (GrB_Index t = 0 ; t < Anz ; t++) Ap[I_tup[t] + 1]++ ;
            for (GrB_Index r = 0 ; r < n_cur ; r++) Ap[r + 1] += Ap[r] ;
            memcpy (cursor, Ap, n_cur * sizeof (GrB_Index)) ;
            for (GrB_Index t = 0 ; t < Anz ; t++)
            {
                GrB_Index r = I_tup[t] ;
                GrB_Index dst = cursor[r]++ ;
                Aj[dst] = J_tup[t] ;
                Ax[dst] = X_tup[t] ;
            }
        }
#endif

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

                GrB_Index row_begin = Ap[i] ;
                GrB_Index row_end   = Ap[i + 1] ;
                if (row_begin == row_end) continue ;

                // Temporarily remove i from community ci.
                k_comm[ci] -= ki ;

                GrB_Index ndirty = 0 ;
                for (GrB_Index t = row_begin ; t < row_end ; t++)
                {
                    GrB_Index j = Aj[t] ;
                    if (j == i) continue ;          // skip self-loop (in A_agg)
                    int64_t cj = c_arr[j] ;
                    if (!dirty[cj])
                    {
                        dirty[cj]            = 1 ;
                        dirty_list[ndirty++] = (GrB_Index) cj ;
                        T_local[cj]          = 0.0 ;
                    }
                    T_local[cj] += Ax[t] ;
                }

                double  T_ci      = dirty[ci] ? T_local[ci] : 0.0 ;
                double  score_ci  = T_ci - ki * k_comm[ci] / two_m ;
                double  best_score = score_ci ;
                int64_t best_c     = ci ;

                for (GrB_Index d = 0 ; d < ndirty ; d++)
                {
                    int64_t c_cand = (int64_t) dirty_list[d] ;
                    if (c_cand == ci) continue ;
                    double score = T_local[c_cand] - ki * k_comm[c_cand] / two_m ;
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

                GrB_Index row_begin = Ap[i] ;
                GrB_Index row_end   = Ap[i + 1] ;
                if (row_begin == row_end) continue ;

                k_ref_comm[ci_ref] -= ki ;

                GrB_Index ndirty = 0 ;
                for (GrB_Index t = row_begin ; t < row_end ; t++)
                {
                    GrB_Index j = Aj[t] ;
                    if (j == i) continue ;              // skip self-loop
                    if (c_p1[j] != pi) continue ;       // cross-parent: skip

                    int64_t cj_ref = c_ref[j] ;
                    if (!dirty[cj_ref])
                    {
                        dirty[cj_ref]        = 1 ;
                        dirty_list[ndirty++] = (GrB_Index) cj_ref ;
                        T_local[cj_ref]      = 0.0 ;
                    }
                    T_local[cj_ref] += Ax[t] ;
                }

                double  T_ci_ref    = dirty[ci_ref] ? T_local[ci_ref] : 0.0 ;
                double  score_ci_ref = T_ci_ref - ki * k_ref_comm[ci_ref] / two_m ;
                double  best_score  = score_ci_ref ;
                int64_t best_c_ref  = ci_ref ;

                for (GrB_Index d = 0 ; d < ndirty ; d++)
                {
                    int64_t c_cand = (int64_t) dirty_list[d] ;
                    if (c_cand == ci_ref) continue ;
                    double score = T_local[c_cand] - ki * k_ref_comm[c_cand] / two_m ;
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
            o_comm[i] = c_ref[o_comm[i]] ;
        }

        //----------------------------------------------------------------------
        // Reload A_cur from the container before any further GraphBLAS use
        // (Phase 3 mxm or next-level unload).  Ownership of Ap/Aj/Ax returns
        // to GraphBLAS; we null our pointers so LG_FREE_WORK won't double-free.
        // No-op on the v9 fallback (A_cur was never unloaded).
        //----------------------------------------------------------------------

#if LAGR_LEIDEN_USE_CONTAINER
        GRB_TRY (GxB_Vector_load (cont->p, (void **) &Ap, pty,
            pn, psz, ph, NULL)) ;
        Ap = NULL ;
        GRB_TRY (GxB_Vector_load (cont->i, (void **) &Aj, ity,
            in_, isz, ih, NULL)) ;
        Aj = NULL ;
        GRB_TRY (GxB_Vector_load (cont->x, (void **) &Ax, xty,
            xn, xsz, xh, NULL)) ;
        Ax = NULL ;
        GRB_TRY (GxB_load_Matrix_from_Container (A_cur, cont, NULL)) ;
#endif

        //----------------------------------------------------------------------
        // PHASE 3: Aggregation — build coarsened graph if communities merged
        //----------------------------------------------------------------------

        if (K_ref < n_cur)
        {
            outer_changed = true ;

            // S_mat: n_cur × K_ref indicator matrix; S[i, c_ref[i]] = 1 for
            // every i.  All values are 1.0, so use GxB_Matrix_build_Scalar
            // and a shared scalar instead of materializing a 1.0-array.
            //   rows: iota (precomputed [0..n-1], reused)
            //   cols: c_ref reinterpret-cast to GrB_Index*.  c_ref values
            //         are non-negative community labels in [0, K_ref); on
            //         all targeted platforms int64_t and uint64_t share the
            //         same width and representation for non-negative values.
            GRB_TRY (GrB_Matrix_new (&S_mat, GrB_FP64, n_cur, K_ref)) ;
            GRB_TRY (GxB_Matrix_build_Scalar (S_mat, iota,
                (GrB_Index *) c_ref, one_scalar, n_cur)) ;

            // A_temp = A_cur * S  (n_cur × K_ref)
            GRB_TRY (GrB_Matrix_new (&A_temp, GrB_FP64, n_cur, K_ref)) ;
            GRB_TRY (GrB_mxm (A_temp, NULL, NULL,
                GrB_PLUS_TIMES_SEMIRING_FP64, A_cur, S_mat, NULL)) ;

            // A_new = S^T * A_temp  (K_ref × K_ref)
            GRB_TRY (GrB_Matrix_new (&A_new, GrB_FP64, K_ref, K_ref)) ;
            GRB_TRY (GrB_mxm (A_new, NULL, NULL,
                GrB_PLUS_TIMES_SEMIRING_FP64, S_mat, A_temp, GrB_DESC_T0)) ;

            GrB_free (&S_mat) ;  S_mat  = NULL ;
            GrB_free (&A_temp) ; A_temp = NULL ;
            GrB_free (&A_agg) ;  // free previous level's aggregate graph
            A_agg  = A_new ;
            A_new  = NULL ;     // ownership transferred to A_agg
            A_cur  = A_agg ;
            n_cur  = K_ref ;
        }
        // K_ref == n_cur: no communities merged this level → converged.
    }

    //--------------------------------------------------------------------------
    // Build output GrB_Vector from o_comm with move semantics: hand the
    // o_comm buffer directly to GraphBLAS (no copy) and null our pointer so
    // LG_FREE_WORK doesn't double-free.  o_comm values are already relabeled
    // 0..K_final-1 from the last iteration; the loaded vector is "full"
    // (every index has a value), so set sparsity hint accordingly.
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (c_handle, GrB_INT64, n)) ;
#if LAGR_LEIDEN_USE_CONTAINER
    GRB_TRY (GrB_set (*c_handle, GxB_FULL, GxB_SPARSITY_CONTROL)) ;
    GRB_TRY (GxB_Vector_load (*c_handle, (void **) &o_comm, GrB_INT64,
        n, n * sizeof (int64_t), GrB_DEFAULT, NULL)) ;
    o_comm = NULL ;     // ownership transferred to *c_handle
#else
    GRB_TRY (GrB_Vector_build_INT64 (*c_handle, iota, o_comm, n,
        GrB_FIRST_INT64)) ;
#endif

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
