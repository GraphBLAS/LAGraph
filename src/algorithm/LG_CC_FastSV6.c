//------------------------------------------------------------------------------
// LG_CC_FastSV6: connected components
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Code is based on the algorithm described in the following paper:
// Zhang, Azad, Hu. FastSV: A Distributed-Memory Connected Component
// Algorithm with Fast Convergence (SIAM PP20)

// A subsequent update to the algorithm is here (which might not be reflected
// in this code):
// Yongzhe Zhang, Ariful Azad, Aydin Buluc: Parallel algorithms for finding
// connected components using linear algebra. J. Parallel Distributed Comput.
// 144: 14-27 (2020).

// Modified by Tim Davis, Texas A&M University: revised Reduce_assign to use
// purely GrB* and GxB* methods and the matrix C.  Added warmup phase.  Changed
// to use GxB pack/unpack instead of GxB import/export.  Converted to use the
// LAGraph_Graph object.  Exploiting iso property for the temporary matrices
// C and T.

// The input graph G must be undirected, or directed and with an adjacency
// matrix that has a symmetric structure.  Self-edges (diagonal entries) are
// OK, and are ignored.  The values and type of A are ignored; just its
// structure is accessed.

// TODO: This function must not be called by multiple user threads at the same
// time on the same graph G, since it unpacks G->A and then packs it back when
// done.  G->A is unchanged when the function returns, but during execution
// G->A is empty.  This will be fixed once the TODOs are finished below, and
// G->A will then become a truly read-only object (assuming GrB_wait (G->A)
// has been done first).

#define LAGraph_FREE_ALL ;
#include "LG_internal.h"

#if !LG_VANILLA
#if (! LG_SUITESPARSE )
#error "SuiteSparse:GraphBLAS v6.0.1 or later required"
#endif

//==============================================================================
// fastsv: find the components of a graph
//==============================================================================

static inline GrB_Info fastsv
(
    GrB_Matrix A,           // adjacency matrix, G->A or a subset of G->A
    GrB_Vector parent,      // parent vector
    GrB_Vector mngp,        // min neighbor grandparent
    GrB_Vector *gp,         // grandparent
    GrB_Vector *gp_new,     // new grandparent (swapped with gp)
    GrB_Vector t,           // workspace
    GrB_BinaryOp eq,        // GrB_EQ_(integer type)
    GrB_BinaryOp min,       // GrB_MIN_(integer type)
    GrB_Semiring min_2nd,   // GrB_MIN_SECOND_(integer type)
    GrB_Matrix C,           // C(i,j) present if i = Px (j)
    GrB_Index **Cp,         // 0:n, size n+1
    GrB_Index **Px,         // Px: non-opaque copy of parent vector, size n
    void **Cx,              // size 1, contents not accessed
    char *msg
)
{
    GrB_Index n ;
    GrB_TRY (GrB_Vector_size (&n, parent)) ;
    GrB_Index Cp_size = (n+1) * sizeof (GrB_Index) ;
    GrB_Index Ci_size = n * sizeof (GrB_Index) ;
    GrB_Index Cx_size = sizeof (bool) ;
    bool iso = true, jumbled = false, done = false ;

    while (true)
    {

        //----------------------------------------------------------------------
        // hooking & shortcutting
        //----------------------------------------------------------------------

        // mngp = min (mngp, A*gp) using the MIN_SECOND semiring
        GrB_TRY (GrB_mxv (mngp, NULL, min, min_2nd, A, *gp, NULL)) ;

        //----------------------------------------------------------------------
        // parent = min (parent, C*mngp) where C(i,j) is present if i=Px(j)
        //----------------------------------------------------------------------

        // Reduce_assign: The Px array of size n is the non-opaque copy of the
        // parent vector, where i = Px [j] if the parent of node j is node i.
        // It can thus have duplicates.  The vectors parent and mngp are full
        // (all entries present).  This function computes the following, which
        // is done explicitly in the Reduce_assign function in LG_CC_Boruvka:
        //
        //      for (j = 0 ; j < n ; j++)
        //      {
        //          uint64_t i = Px [j] ;
        //          parent [i] = min (parent [i], mngp [j]) ;
        //      }
        //
        // If C(i,j) is present where i == Px [j], then this can be written as:
        //
        //      parent = min (parent, C*mngp)
        //
        // when using the min_2nd semiring.  This can be done efficiently
        // because C can be constructed in O(1) time and O(1) additional space
        // (not counting the prior Cp, Px, and Cx arrays), when using the
        // SuiteSparse pack/unpack move constructors.  The min_2nd semiring
        // ignores the values of C and operates only on the structure, so its
        // values are not relevant.  Cx is thus chosen as a GrB_BOOL array of
        // size 1 where Cx [0] = false, so the all entries present in C are
        // equal to false.

        // pack Cp, Px, and Cx into a matrix C with C(i,j) present if Px(j) == i
        GrB_TRY (GxB_Matrix_pack_CSC (C, Cp, /* Px is Ci: */ Px, Cx,
            Cp_size, Ci_size, Cx_size, iso, jumbled, NULL)) ;

        // parent = min (parent, C*mngp) using the MIN_SECOND semiring
        GrB_TRY (GrB_mxv (parent, NULL, min, min_2nd, C, mngp, NULL)) ;

        // unpack the contents of C, to make Px available to this method again.
        GrB_TRY (GxB_Matrix_unpack_CSC (C, Cp, Px, Cx,
            &Cp_size, &Ci_size, &Cx_size, &iso, &jumbled, NULL)) ;

        //----------------------------------------------------------------------
        // parent = min (parent, mngp, gp)
        //----------------------------------------------------------------------

        GrB_TRY (GrB_eWiseAdd (parent, NULL, min, min, mngp, *gp, NULL)) ;

        //----------------------------------------------------------------------
        // calculate grandparent: gp_new = parent (parent), and extract Px
        //----------------------------------------------------------------------

        // if parent is uint32, GraphBLAS typecasts to uint64 for Px.
        GrB_TRY (GrB_Vector_extractTuples (NULL, *Px, &n, parent)) ;
        GrB_TRY (GrB_extract (*gp_new, NULL, NULL, parent, *Px, n, NULL)) ;

        //----------------------------------------------------------------------
        // terminate if gp and gp_new are the same
        //----------------------------------------------------------------------

        GrB_TRY (GrB_eWiseMult (t, NULL, NULL, eq, *gp_new, *gp, NULL)) ;
        GrB_TRY (GrB_reduce (&done, NULL, GrB_LAND_MONOID_BOOL, t, NULL)) ;
        if (done) break ;

        // swap gp and gp_new
        GrB_Vector s = (*gp) ; (*gp) = (*gp_new) ; (*gp_new) = s ;
    }
    return (GrB_SUCCESS) ;
}

//==============================================================================
// LG_CC_FastSV6
//==============================================================================

// The output of LG_CC_FastSV* is a vector component, where component(i)=r if
// node i is in the connected compononent whose representative is node r.  If r
// is a representative, then component(r)=r.  The number of connected
// components in the graph G is the number of representatives.

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK                   \
{                                           \
    LAGraph_Free ((void **) &Tp) ;          \
    LAGraph_Free ((void **) &Tj) ;          \
    LAGraph_Free ((void **) &Tx) ;          \
    LAGraph_Free ((void **) &Cp) ;          \
    LAGraph_Free ((void **) &Px) ;          \
    LAGraph_Free ((void **) &Cx) ;          \
    LAGraph_Free ((void **) &ht_key) ;      \
    LAGraph_Free ((void **) &ht_count) ;    \
    LAGraph_Free ((void **) &count) ;       \
    LAGraph_Free ((void **) &range) ;       \
    GrB_free (&C) ;                         \
    GrB_free (&T) ;                         \
    GrB_free (&t) ;                         \
    GrB_free (&y) ;                         \
    GrB_free (&gp) ;                        \
    GrB_free (&mngp) ;                      \
    GrB_free (&gp_new) ;                    \
}

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL                    \
{                                           \
    LAGraph_FREE_WORK ;                     \
    GrB_free (&parent) ;                    \
}

#endif

int LG_CC_FastSV6           // SuiteSparse:GraphBLAS method, with GxB extensions
(
    // output
    GrB_Vector *component,  // component(i)=r if node is in the component r
    // inputs
    LAGraph_Graph G,        // input graph
    char *msg
)
{

#if LG_VANILLA
    LG_CHECK (0, GrB_NOT_IMPLEMENTED, "SuiteSparse required for this method") ;
#else

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    int64_t *range = NULL ;
    GrB_Index n, nvals, Cp_size = 0, *ht_key = NULL, *Px = NULL, *Cp = NULL,
        *count = NULL, *Tp = NULL, *Tj = NULL ;
    GrB_Vector parent = NULL, gp_new = NULL, mngp = NULL, gp = NULL, t = NULL,
        y = NULL ;
    GrB_Matrix T = NULL, C = NULL ;
    void *Tx = NULL, *Cx = NULL ;
    int *ht_count = NULL ;

    LG_CHECK (LAGraph_CheckGraph (G, msg), GrB_INVALID_OBJECT,
        "graph is invalid") ;
    LG_CHECK (component == NULL, GrB_NULL_POINTER, "component is NULL") ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_structure_is_symmetric == LAGRAPH_TRUE))
    {
        // A must be symmetric
        ;
    }
    else
    {
        // A must not be unsymmetric
        LG_CHECK (false, GrB_INVALID_VALUE, "input must be symmetric") ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;

    // determine the integer type, operators, and semirings to use 
    GrB_Type Uint, Int ;
    GrB_IndexUnaryOp ramp ;
    GrB_Semiring min_2nd, min_2ndi ;
    GrB_BinaryOp min, eq, imin ;
    #ifdef COVERAGE
    // Just for test coverage, use 64-bit ints for n > 100.  Do not use this
    // rule in production!
    #define NBIG 100
    #else
    // For production use: 64-bit integers if n > 2^31
    #define NBIG INT32_MAX
    #endif
    if (n > NBIG)
    {
        // use 64-bit integers throughout
        Uint = GrB_UINT64 ;
        Int  = GrB_INT64  ;
        ramp = GrB_ROWINDEX_INT64 ;
        min  = GrB_MIN_UINT64 ;
        imin = GrB_MIN_INT64 ;
        eq   = GrB_EQ_UINT64 ;
        min_2nd  = GrB_MIN_SECOND_SEMIRING_UINT64 ;
        min_2ndi = GxB_MIN_SECONDI_INT64 ;
    }
    else
    {
        // use 32-bit integers, except for Px and for constructing the matrix C
        Uint = GrB_UINT32 ;
        Int  = GrB_INT32  ;
        ramp = GrB_ROWINDEX_INT32 ;
        min  = GrB_MIN_UINT32 ;
        imin = GrB_MIN_INT32 ;
        eq   = GrB_EQ_UINT32 ;
        min_2nd  = GrB_MIN_SECOND_SEMIRING_UINT32 ;
        min_2ndi = GxB_MIN_SECONDI_INT32 ;
    }

    // FASTSV_SAMPLES: number of samples to take from each row A(i,:).
    // Sampling is used if the average degree is > 8 and if n > 1024.
    #define FASTSV_SAMPLES 4
    bool sampling = (nvals > n * FASTSV_SAMPLES * 2 && n > 1024) ;

// [ TODO: nthreads will not be needed once GxB_select with a GxB_RankUnaryOp
// and a new GxB_extract are added to SuiteSparse:GraphBLAS.
    // determine # of threads to use
    int nthreads ;
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads, NULL)) ;
    nthreads = LAGraph_MIN (nthreads, n / 16) ;
    nthreads = LAGraph_MAX (nthreads, 1) ;
// ]

    Cx = (void *) LAGraph_Calloc (1, sizeof (bool)) ;
    Px = (GrB_Index *) LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    LG_CHECK (Px == NULL || Cx == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;

    // create Cp = 0:n (always 64-bit) and the empty C matrix
    GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;
    GrB_TRY (GrB_Vector_new (&t, GrB_INT64, n+1)) ;
    GrB_TRY (GrB_assign (t, NULL, NULL, 0, GrB_ALL, n+1, NULL)) ;
    GrB_TRY (GrB_apply (t, NULL, NULL, GrB_ROWINDEX_INT64, t, 0, NULL)) ;
    GrB_TRY (GxB_Vector_unpack_Full (t, (void **) &Cp, &Cp_size, NULL, NULL)) ;
    GrB_TRY (GrB_free (&t)) ;

    //--------------------------------------------------------------------------
    // warmup: parent = min (0:n-1, A*1) using the MIN_SECONDI semiring
    //--------------------------------------------------------------------------

    // y (i) = min (i, j) for all entries A(i,j).  This warmup phase takes only
    // O(n) time, because of how the MIN_SECONDI semiring is implemented in
    // SuiteSparse:GraphBLAS.  A is held by row, and the first entry in A(i,:)
    // is the minimum index j, so only the first entry in A(i,:) needs to be
    // considered for each row i.

    GrB_TRY (GrB_Vector_new (&t, Int, n)) ;
    GrB_TRY (GrB_Vector_new (&y, Int, n)) ;
    GrB_TRY (GrB_assign (t, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_assign (y, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_apply (y, NULL, NULL, ramp, y, 0, NULL)) ;
    GrB_TRY (GrB_mxv (y, NULL, imin, min_2ndi, A, t, NULL)) ;
    GrB_TRY (GrB_free (&t)) ;

    // The typecast from Int to Uint is required because the ROWINDEX operator
    // and MIN_SECONDI do not work in the UINT* domains, as built-in operators.
    // parent = (Uint) y
    GrB_TRY (GrB_Vector_new (&parent, Uint, n)) ;
    GrB_TRY (GrB_assign (parent, NULL, NULL, y, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_free (&y)) ;

    // copy parent into gp, mngp, and Px.  Px is a non-opaque 64-bit copy of the
    // parent GrB_Vector.  The Px array is always of type GrB_Index since it
    // must be used as the input array for extractTuples and as Ci for pack_CSR.
    // If parent is uint32, GraphBLAS typecasts it to the uint64 Px array.
    GrB_TRY (GrB_Vector_extractTuples (NULL, Px, &n, parent)) ;
    GrB_TRY (GrB_Vector_dup (&gp, parent)) ;
    GrB_TRY (GrB_Vector_dup (&mngp, parent)) ;
    GrB_TRY (GrB_Vector_new (&gp_new, Uint, n)) ;
    GrB_TRY (GrB_Vector_new (&t, GrB_BOOL, n)) ;

    //--------------------------------------------------------------------------
    // sample phase
    //--------------------------------------------------------------------------

    if (sampling)
    {

// [ TODO: GxB_select, using a new operator: GxB_RankUnaryOp, will do all this,
// with GxB_Matrix_select_RankOp_Scalar with operator GxB_LEASTRANK and a
// GrB_Scalar input equal to FASTSV_SAMPLES.  Built-in operators will be,
// (where y is INT64):
//
//      GxB_LEASTRANK (aij, i, j, k, d, y): select if aij has rank k <= y
//      GxB_NLEASTRANK: select if aij has rank k > y
//      GxB_GREATESTRANK (...) select if aij has rank k >= (d-y) where
//          d = # of entries in A(i,:).
//      GxB_NGREATESTRANK (...): select if aij has rank k < (d-y)
// and perhaps other operators such as:
//      GxB_LEASTRELRANK (...): select aij if rank k <= y*d where y is double
//      GxB_GREATESTRELRANK (...): select aij rank k > y*d where y is double
//
// By default, the rank of aij is its relative position as the kth entry in its
// row (from "left" to "right").  If a new GxB setting in the descriptor is
// set, then k is the relative position of aij as the kth entry in its column.
// The default would be that the rank is the position of aij in its row A(i,:).

// Other:
//      give me 3 random items from the row (y = 3)
//      give me the 4 biggest *values* in each row (y = 4)

// mxv:
//      C = A*diag(D)

        //----------------------------------------------------------------------
        // unpack A in CSR format
        //----------------------------------------------------------------------

        void *Ax ;
        GrB_Index *Ap, *Aj, Ap_size, Aj_size, Ax_size ;
        bool A_jumbled, A_iso ;
        GrB_TRY (GxB_Matrix_unpack_CSR (A, &Ap, &Aj, &Ax,
            &Ap_size, &Aj_size, &Ax_size, &A_iso, &A_jumbled, NULL)) ;

        //----------------------------------------------------------------------
        // allocate workspace, including space to construct T
        //----------------------------------------------------------------------

        GrB_Index Tp_size = (n+1) * sizeof (GrB_Index) ;
        GrB_Index Tj_size = nvals * sizeof (GrB_Index) ;
        GrB_Index Tx_size = sizeof (bool) ;
        Tp = (GrB_Index *) LAGraph_Malloc (n+1, sizeof (GrB_Index)) ;
        Tj = (GrB_Index *) LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;
        Tx = (bool *) LAGraph_Calloc (1, sizeof (bool)) ;
        range = (int64_t *) LAGraph_Malloc (nthreads + 1, sizeof (int64_t)) ;
        count = (GrB_Index *) LAGraph_Calloc (nthreads + 1, sizeof (GrB_Index));
        LG_CHECK (Tp == NULL || Tj == NULL || Tx == NULL || range == NULL
            || count == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;

        //----------------------------------------------------------------------
        // define parallel tasks to construct T
        //----------------------------------------------------------------------

        // thread tid works on rows range[tid]:range[tid+1]-1 of A and T
        for (int tid = 0 ; tid <= nthreads ; tid++)
        {
            range [tid] = (n * tid + nthreads - 1) / nthreads ;
        }

        //----------------------------------------------------------------------
        // determine the number entries to be constructed in T for each thread
        //----------------------------------------------------------------------

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                int64_t deg = Ap [i + 1] - Ap [i] ;
                count [tid + 1] += LAGraph_MIN (FASTSV_SAMPLES, deg) ;
            }
        }

        //----------------------------------------------------------------------
        // count = cumsum (count)
        //----------------------------------------------------------------------

        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            count [tid + 1] += count [tid] ;
        }

        //----------------------------------------------------------------------
        // construct T
        //----------------------------------------------------------------------

        // T (i,:) consists of the first FASTSV_SAMPLES of A (i,:).

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = count [tid] ;
            Tp [range [tid]] = p ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                // construct T (i,:) from the first entries in A (i,:)
                for (int64_t j = 0 ;
                    j < FASTSV_SAMPLES && Ap [i] + j < Ap [i + 1] ; j++)
                {
                    Tj [p++] = Aj [Ap [i] + j] ;
                }
                Tp [i + 1] = p ;
            }
        }

        //----------------------------------------------------------------------
        // import the result into the GrB_Matrix T
        //----------------------------------------------------------------------

        GrB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;
        GrB_TRY (GxB_Matrix_pack_CSR (T, &Tp, &Tj, &Tx, Tp_size, Tj_size,
            Tx_size, /* T is iso: */ true, A_jumbled, NULL)) ;

// ] TODO: the above will all be done as a single call to GxB_select.

        //----------------------------------------------------------------------
        // find the connected components of T
        //----------------------------------------------------------------------

        GrB_TRY (fastsv (T, parent, mngp, &gp, &gp_new, t, eq, min, min_2nd,
            C, &Cp, &Px, &Cx, msg)) ;

        //----------------------------------------------------------------------
        // use sampling to estimate the largest connected component in T
        //----------------------------------------------------------------------

        // The sampling below computes an estimate of the mode of the parent
        // vector, the contents of which are currently in the non-opaque Px
        // array.

        // FIXME: brittle:  864, 1024

        // hash table size must be a power of 2
        #define HASH_SIZE 1024
        // number of samples to insert into the hash table
        #define HASH_SAMPLES 864
        #define HASH(x) (((x << 4) + x) & (HASH_SIZE-1))
        #define NEXT(x) ((x + 23) & (HASH_SIZE-1))

        // allocate and initialize the hash table
        ht_key = (GrB_Index *) LAGraph_Malloc (HASH_SIZE, sizeof (GrB_Index)) ;
        ht_count = (int *) LAGraph_Calloc (HASH_SIZE, sizeof (int)) ;
        LG_CHECK (ht_key == NULL || ht_count == NULL, GrB_OUT_OF_MEMORY,
            "out of memory") ;
        for (int k = 0 ; k < HASH_SIZE ; k++)
        {
            ht_key [k] = UINT64_MAX ;
        }

        // hash the samples and find the most frequent entry
        uint64_t seed = n ;         // random number seed
        int64_t key = -1 ;          // most frequent entry
        int max_count = 0 ;         // frequency of most frequent entry
        for (int64_t k = 0 ; k < HASH_SAMPLES ; k++)
        {
            // select an entry from Px at random
            GrB_Index x = Px [LAGraph_Random60 (&seed) % n] ;
            // find x in the hash table
            GrB_Index h = HASH (x) ;
            while (ht_key [h] != UINT64_MAX && ht_key [h] != x) h = NEXT (h) ;
            // add x to the hash table
            ht_key [h] = x ;
            ht_count [h]++ ;
            // keep track of the most frequent value
            if (ht_count [h] > max_count)
            {
                key = ht_key [h] ;
                max_count = ht_count [h] ;
            }
        }

        //----------------------------------------------------------------------
        // compact the largest connected component in A
        //----------------------------------------------------------------------

        // Construct a new matrix T from the input matrix A (the matrix A is
        // not changed). The key node is the representative of the (estimated)
        // largest component.  T is constructed as a copy of A, except:
        // (1) all edges A(i,:) for nodes i in the key component deleted, and
        // (2) for nodes i not in the key component, A(i,j) is deleted if 
        //     j is in the key component.
        // (3) If A(i,:) has any deletions from (2), T(i,key) is added to T.

// [ TODO: replace this with GxB_extract with GrB_Vector index arrays.
// See https://github.com/GraphBLAS/graphblas-api-c/issues/67 .
// This method will not insert the new entries T(i,key) for rows i that have
// had entries deleted.  That can be done with GrB_assign, with an n-by-1 mask
// M computed from the before-and-after row degrees of A and T:
// M = (parent != key) && (row_degree(T) < row_degree(A))
// J [0] = key.
// GxB_Matrix_subassign_BOOL (T, M, NULL, true, GrB_ALL, n, J, 1, NULL)
// or with
// GrB_Col_assign (T, M, NULL, t, GrB_ALL, j, NULL) with an all-true
// vector t.

        // unpack T to reuse the space (all content is overwritten below)
        bool T_jumbled, T_iso ;
        GrB_TRY (GxB_Matrix_unpack_CSR (T, &Tp, &Tj, &Tx, &Tp_size, &Tj_size,
            &Tx_size, &T_iso, &T_jumbled, NULL)) ;

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = Ap [range [tid]] ;
            // thread tid scans A (range [tid]:range [tid+1]-1,:),
            // and constructs T(i,:) for all rows in this range.
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                int64_t pi = Px [i] ;   // pi = parent (i)
                Tp [i] = p ;            // start the construction of T(i,:)
                // T(i,:) is empty if pi == key
                if (pi != key)
                {
                    // scan A(i,:)
                    for (GrB_Index pS = Ap [i] ; pS < Ap [i+1] ; pS++)
                    {
                        // get A(i,j)
                        int64_t j = Aj [pS] ;
                        if (Px [j] != key)
                        {
                            // add the entry T(i,j) to T, but skip it if
                            // Px [j] is equal to key
                            Tj [p++] = j ;
                        }
                    }
                    // Add the entry T(i,key) if there is room for it in T(i,:);
                    // if and only if node i is adjacent to a node j in the
                    // largest component.  The only way there can be space if
                    // at least one T(i,j) appears with Px [j] equal to the key
                    // (that is, node j is in the largest connected component,
                    // key == Px [j].  One of these j's can then be replaced
                    // with the key.  If node i is not adjacent to any node in
                    // the largest component, then there is no space in T(i,:)
                    // and no new edge to the largest component is added.
                    if (p - Tp [i] < Ap [i+1] - Ap [i])
                    {
                        Tj [p++] = key ;
                    }
                }
            }
            // count the number of entries inserted into T by this thread
            count [tid] = p - Tp [range [tid]] ;
        }

        // Compact empty space out of Tj not filled in from the above phase.
        nvals = 0 ;
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            memcpy (Tj + nvals, Tj + Tp [range [tid]],
                sizeof (GrB_Index) * count [tid]) ;
            nvals += count [tid] ;
            count [tid] = nvals - count [tid] ;
        }

        // Compact empty space out of Tp
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = Tp [range [tid]] ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                Tp [i] -= p - count [tid] ;
            }
        }

        // finalize T
        Tp [n] = nvals ;

        // pack T for the final phase
        GrB_TRY (GxB_Matrix_pack_CSR (T, &Tp, &Tj, &Tx, Tp_size, Tj_size,
            Tx_size, T_iso, /* T is now jumbled */ true, NULL)) ;

        // pack A (unchanged since last unpack); this is the original G->A.
        GrB_TRY (GxB_Matrix_pack_CSR (A, &Ap, &Aj, &Ax, Ap_size, Aj_size,
            Ax_size, A_iso, A_jumbled, NULL)) ;

// ].  The unpack/pack of A into Ap, Aj, Ax will not be needed, and G->A
// will become truly a read-only matrix.

        // final phase uses the pruned matrix T
        A = T ;
    }

    //--------------------------------------------------------------------------
    // check for quick return
    //--------------------------------------------------------------------------

    // The sample phase may have already found that G->A has a single component,
    // in which case the matrix A is now empty.

    if (nvals == 0)
    {
        (*component) = parent ;
        LAGraph_FREE_WORK ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // final phase
    //--------------------------------------------------------------------------

    GrB_TRY (fastsv (A, parent, mngp, &gp, &gp_new, t, eq, min, min_2nd,
        C, &Cp, &Px, &Cx, msg)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*component) = parent ;
    LAGraph_FREE_WORK ;
    return (GrB_SUCCESS) ;
#endif
}
