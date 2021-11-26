//------------------------------------------------------------------------------
// LG_CC_FastSV6: connected components
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Code is based on the algorithm described in the following paper
// Zhang, Azad, Hu. FastSV: FastSV: A Distributed-Memory Connected Component
// Algorithm with Fast Convergence (SIAM PP20)

// A subsequent update to the algorithm is here (which might not be reflected
// in this code):
// Yongzhe Zhang, Ariful Azad, Aydin Buluc: Parallel algorithms for finding
// connected components using linear algebra. J. Parallel Distributed Comput.
// 144: 14-27 (2020).

// Modified by Tim Davis, Texas A&M University

// The input graph G must be undirected, or directed and with an adjacency
// matrix that has a symmetric structure.  Self-edges (diagonal entries) are
// OK, and are ignored.  The values and type of A are ignored; just its
// structure is accessed.

// This function should not be called by multiple user threads on the same
// graph G, since it unpacks G->A and then packs it back when done.  G->A is
// unchanged when the function returns, but during execution G->A is empty.

#define LAGraph_FREE_ALL ;
#include "LG_internal.h"

#if !LG_VANILLA
#if (! LG_SUITESPARSE )
#error "SuiteSparse:GraphBLAS v6.0.0 or later required"
#endif

//------------------------------------------------------------------------------
// Reduce_assign:  w (Px) += s, using MIN as the "+=" accum operator
//------------------------------------------------------------------------------

// The Px array of size n is the non-opaque version of the parent vector, where
// i = Px [j] if the parent of node j is node i.  It can thus have duplicates.
// The vectors w and s are full (all entries present).  This function computes
// the following, which is done explicitly in the Reduce_assign function in
// LG_CC_Boruvka:
//
//      for (j = 0 ; j < n ; j++)
//      {
//          uint64_t i = Px [j] ;
//          w [i] = min (w [i], s [j]) ;
//      }
//
// If C(i,j) = true where i == Px [j], then this can be written as:
//
//      w = min (w, C*s)
//
// when using the min_second semiring.  This can be done efficiently where
// because C can be constructed in O(1) time and O(1) additional space (not
// counting the prior Cp, Ci, and Cx arrays), when using the SuiteSparse
// pack/unpack move constructors.

static inline GrB_Info Reduce_assign
(
    GrB_Vector w,           // vector of size n, all entries present
    GrB_Vector s,           // vector of size n, all entries present
    GrB_Matrix C,           // boolean matrix of size n-by-n
    GrB_Index **Cp_handle,  // array of size n+1, equal to 0:n
    GrB_Index **Ci_handle,  // Px array of size n, can have duplicates
    bool **Cx_handle,       // array of size 1, equal to true
    char *msg
)
{
    // size of Cp, Ci, and Cx in bytes
    GrB_Index n ;
    GrB_TRY (GrB_Vector_size (&n, w)) ;
    GrB_Index Cp_size = (n+1) * sizeof (GrB_Index) ;
    GrB_Index Ci_size = n * sizeof (GrB_Index) ;
    GrB_Index Cx_size = sizeof (bool) ;

    // pack Cp, Ci, and Cx into a matrix C with C(i,j) = true if Ci(j) == i
    bool iso = true ;
    bool jumbled = false ;
    GrB_TRY (GxB_Matrix_pack_CSC (C, Cp_handle, Ci_handle, (void **) Cx_handle,
        Cp_size, Ci_size, Cx_size, iso, jumbled, NULL)) ;

    // w = min (w, C*s) using the MIN_SECOND semiring
    GrB_TRY (GrB_mxv (w, NULL, GrB_MIN_UINT64,
        GrB_MIN_SECOND_SEMIRING_UINT64, C, s, NULL)) ;

    // unpack the contents of C
    GrB_TRY (GxB_Matrix_unpack_CSC (C, Cp_handle, Ci_handle, (void **)Cx_handle,
        &Cp_size, &Ci_size, &Cx_size, &iso, &jumbled, NULL)) ;

    return (GrB_SUCCESS) ;  // yay! It works!
}

//------------------------------------------------------------------------------
// LG_CC_FastSV6
//------------------------------------------------------------------------------

// The output of LG_CC_FastSV* is a vector component, where
// component(i)=s if node i is in the connected compononent whose
// representative node is node s.  If s is a representative, then
// component(s)=s.  The number of connected components in the graph G is the
// number of representatives.

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK                   \
{                                           \
    LAGraph_Free ((void **) &Tp) ;          \
    LAGraph_Free ((void **) &Tj) ;          \
    LAGraph_Free ((void **) &Tx) ;          \
    LAGraph_Free ((void **) &Cp) ;          \
    LAGraph_Free ((void **) &Cx) ;          \
    LAGraph_Free ((void **) &Px) ;          \
    LAGraph_Free ((void **) &ht_key) ;      \
    LAGraph_Free ((void **) &ht_count) ;    \
    LAGraph_Free ((void **) &count) ;       \
    LAGraph_Free ((void **) &range) ;       \
    GrB_free (&T) ;                         \
    GrB_free (&t) ;                         \
    GrB_free (&y) ;                         \
    GrB_free (&gp) ;                        \
    GrB_free (&mngp) ;                      \
    GrB_free (&gp_new) ;                    \
    GrB_free (&mod) ;                       \
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
    GrB_Vector *component,  // component(i)=s if node is in the component s
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

    int64_t *ht_key = NULL, *ht_count = NULL, *range = NULL ;
    GrB_Index n, nnz, Cp_size = 0,
        *Px = NULL, *Cp = NULL, *count = NULL, *Tp = NULL, *Tj = NULL ;
    GrB_Vector parent = NULL, gp_new = NULL, mngp = NULL, mod = NULL, gp = NULL,
        t = NULL, y = NULL ;
    GrB_Matrix T = NULL, C = NULL ;
    bool *Cx = NULL ;
    void *Tx = NULL ;

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
    GrB_TRY (GrB_Matrix_nvals (&nnz, A)) ;

    // FASTSV_SAMPLES: number of samples to take from each row A(i,:)
    #define FASTSV_SAMPLES 4
    bool sampling = (n * FASTSV_SAMPLES * 2 < nnz && n > 1024) ;

    // determine # of threads to use
    int nthreads ;
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads, NULL)) ;
    nthreads = LAGraph_MIN (nthreads, n / 16) ;
    nthreads = LAGraph_MAX (nthreads, 1) ;

    GrB_TRY (GrB_Vector_new (&gp_new, GrB_UINT64, n)) ;
    GrB_TRY (GrB_Vector_new (&mod,    GrB_BOOL,   n)) ;

    Cx = (bool *) LAGraph_Malloc (1, sizeof (bool)) ;
    Px = LAGraph_Malloc (n, sizeof (uint64_t)) ;
    LG_CHECK (Px == NULL || Cx == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;
    Cx [0] = true ;

    // create Cp = 0:n and the empty C matrix
    GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;
    GrB_TRY (GrB_Vector_new (&t, GrB_INT64, n+1)) ;
    GrB_TRY (GrB_assign (t, NULL, NULL, 0, GrB_ALL, n+1, NULL)) ;
    GrB_TRY (GrB_apply (t, NULL, NULL, GrB_ROWINDEX_INT64, t, 0, NULL)) ;
    GrB_TRY (GxB_Vector_unpack_Full (t, (void **) &Cp, &Cp_size, NULL, NULL)) ;
    GrB_TRY (GrB_free (&t)) ;

    //--------------------------------------------------------------------------
    // warmup: y = min (0:n-1, A*t) using the MIN_SECONDI semiring
    //--------------------------------------------------------------------------

    // y (i) = min (i, j) for all entries A(i,j).  This warmup phase takes only
    // O(n) time, because of how the MIN_SECONDI semiring is implemented in
    // SuiteSparse:GraphBLAS.  A is held by row, and the first entry in A(i,:)
    // is the minimum index j, so only the first entry in A(i,:) needs to be
    // considered for each row i.

    GrB_TRY (GrB_Vector_new (&t, GrB_INT64, n)) ;
    GrB_TRY (GrB_Vector_new (&y, GrB_INT64, n)) ;
    GrB_TRY (GrB_assign (t, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_assign (y, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_apply (y, NULL, NULL, GrB_ROWINDEX_INT64, y, 0, NULL)) ;
    GrB_TRY (GrB_mxv (y, NULL, GrB_MIN_INT64, GxB_MIN_SECONDI_INT64, A, t,
        NULL)) ;
    GrB_TRY (GrB_free (&t)) ;

    // The typecast is required because the ROWINDEX operator and MIN_SECONDI
    // do not work in the UINT64 domain, as built-in operators.
    // parent = (uint64) y
    GrB_TRY (GrB_Vector_new (&parent, GrB_UINT64, n)) ;
    GrB_TRY (GrB_assign (parent, NULL, NULL, y, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_free (&y)) ;

    // copy parent into gp, mngp, and Px.  Px is a non-opaque copy of the
    // parent GrB_Vector
    GrB_TRY (GrB_Vector_extractTuples (NULL, Px, &n, parent)) ;
    GrB_TRY (GrB_Vector_dup (&gp,   parent)) ;
    GrB_TRY (GrB_Vector_dup (&mngp, parent)) ;

    //--------------------------------------------------------------------------
    // sample phase
    //--------------------------------------------------------------------------

    if (sampling)
    {

        //----------------------------------------------------------------------
        // unpack A in CSR format
        //----------------------------------------------------------------------

        void *Sx ;
        GrB_Index *Sp, *Sj, Sp_size, Sj_size, Sx_size, nvals ;
        bool S_jumbled, S_iso ;
        GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
        GrB_TRY (GxB_Matrix_unpack_CSR (A, &Sp, &Sj, &Sx,
            &Sp_size, &Sj_size, &Sx_size, &S_iso, &S_jumbled, NULL)) ;

        //----------------------------------------------------------------------
        // allocate workspace, including space to construct T
        //----------------------------------------------------------------------

        GrB_Index Tp_size = (n+1) * sizeof (GrB_Index) ;
        GrB_Index Tj_size = nvals * sizeof (GrB_Index) ;
        GrB_Index Tx_size = sizeof (bool) ;
        Tp = LAGraph_Malloc (n+1, sizeof (GrB_Index)) ;
        Tj = LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;
        Tx = LAGraph_Calloc (1, sizeof (bool)) ;
        range = LAGraph_Malloc (nthreads + 1, sizeof (int64_t)) ;
        count = LAGraph_Calloc (nthreads + 1, sizeof (GrB_Index)) ;
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
                int64_t deg = Sp [i + 1] - Sp [i] ;
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
        // TODO: this could be done by GxB_Select, using a new operator.

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = count [tid] ;
            Tp [range [tid]] = p ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                // construct T (i,:) from the first entries in A (i,:)
                for (int64_t j = 0 ;
                    j < FASTSV_SAMPLES && Sp [i] + j < Sp [i + 1] ; j++)
                {
                    Tj [p++] = Sj [Sp [i] + j] ;
                }
                Tp [i + 1] = p ;
            }
        }

        //----------------------------------------------------------------------
        // import the result into the GrB_Matrix T
        //----------------------------------------------------------------------

        GrB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;
        GrB_TRY (GxB_Matrix_pack_CSR (T, &Tp, &Tj, &Tx, Tp_size, Tj_size,
            Tx_size, /* T is iso: */ true, S_jumbled, NULL)) ;

        //----------------------------------------------------------------------
        // find the connected components of T
        //----------------------------------------------------------------------

        bool changing = true ;
        while (changing)
        {
            // hooking & shortcutting
            // mngp = min (mngp, A*gp) using the MIN_SECOND semiring
            GrB_TRY (GrB_mxv (mngp, NULL, GrB_MIN_UINT64,
                GrB_MIN_SECOND_SEMIRING_UINT64, T, gp, NULL)) ;

            // parent = min (parent, C*mngp) where C is C(i,j) = true if i=Px(j)
            GrB_TRY (Reduce_assign (parent, mngp, C, &Cp, &Px, &Cx, msg)) ;

            // parent = min (parent, mngp, gp)
            GrB_TRY (GrB_eWiseAdd (parent, NULL, GrB_MIN_UINT64, GrB_MIN_UINT64,
                mngp, gp, NULL)) ;

            // calculate grandparent: gp_new = parent (parent)
            GrB_TRY (GrB_Vector_extractTuples (NULL, Px, &n, parent)) ;
            GrB_TRY (GrB_extract (gp_new, NULL, NULL, parent, Px, n, NULL)) ;

            // terminate if gp and gp_new are the same
            GrB_TRY (GrB_eWiseMult (mod, NULL, NULL, GrB_NE_UINT64, gp_new,
                gp, NULL)) ;
            GrB_TRY (GrB_reduce (&changing, NULL, GrB_LOR_MONOID_BOOL, mod,
                NULL)) ;

            // swap gp and gp_new
            GrB_Vector t = gp ; gp = gp_new ; gp_new = t ;
        }

        //----------------------------------------------------------------------
        // use samping to estimate the largest connected component in T
        //----------------------------------------------------------------------

        // hash table size must be a power of 2
        #define HASH_SIZE 1024
        // number of samples to insert into the hash table
        #define HASH_SAMPLES 864
        #define HASH(x) (((x << 4) + x) & (HASH_SIZE-1))
        #define NEXT(x) ((x + 23) & (HASH_SIZE-1))

        // allocate and initialize the hash table
        ht_key = LAGraph_Malloc (HASH_SIZE, sizeof (int64_t)) ;
        ht_count = LAGraph_Calloc (HASH_SIZE, sizeof (int64_t)) ;
        LG_CHECK (ht_key == NULL || ht_count == NULL, GrB_OUT_OF_MEMORY,
            "out of memory") ;
        for (int k = 0 ; k < HASH_SIZE ; k++)
        {
            ht_key [k] = -1 ;
        }

        // hash the samples and find the most frequent entry
        uint64_t seed = n ;         // random number seed
        int64_t key = -1 ;          // most frequent entry
        int64_t max_count = 0 ;     // frequency of most frequent entry
        for (int64_t k = 0 ; k < HASH_SAMPLES ; k++)
        {
            // select an entry from Px at random
            int64_t x = Px [LAGraph_Random60 (&seed) % n] ;
            // find x in the hash table
            int64_t h = HASH (x) ;
            while (ht_key [h] != -1 && ht_key [h] != x)
            {
                h = NEXT (h) ;
            }
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
        // compact the largest connected component in T
        //----------------------------------------------------------------------

        // TODO: replace this with GrB_select and GrB_assign

        // unpack T to resuse the space (all content is overwritten below)
        bool T_jumbled, T_iso ;
        GrB_TRY (GxB_Matrix_unpack_CSR (T, &Tp, &Tj, &Tx, &Tp_size, &Tj_size,
            &Tx_size, &T_iso, &T_jumbled, NULL)) ;

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = Sp [range [tid]] ;
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
                    for (GrB_Index pS = Sp [i] ; pS < Sp [i+1] ; pS++)
                    {
                        // get A(i,j)
                        int64_t j = Sj [pS] ;
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
                    // and no new edge to the largest compenent is added.
                    if (p - Tp [i] < Sp [i+1] - Sp [i])
                    {
                        Tj [p++] = key ;
                    }
                }
            }
            // count the number of entries inserted into T by this thread?
            count [tid] = p - Tp [range [tid]] ;
        }

        // Compact empty space out of Tj not filled in from the above phase.
        nnz = 0 ;
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            memcpy (Tj + nnz, Tj + Tp [range [tid]],
                sizeof (GrB_Index) * count [tid]) ;
            nnz += count [tid] ;
            count [tid] = nnz - count [tid] ;
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
        Tp [n] = nnz ;

        // pack T for the final phase
        GrB_TRY (GxB_Matrix_pack_CSR (T, &Tp, &Tj, &Tx, Tp_size, Tj_size,
            Tx_size, T_iso, /* T is now jumbled */ true, NULL)) ;

        // pack A (unchanged since last unpack)
        GrB_TRY (GxB_Matrix_pack_CSR (A, &Sp, &Sj, &Sx, Sp_size, Sj_size,
            Sx_size, S_iso, S_jumbled, NULL)) ;

        // final phase uses the pruned matrix T
        A = T ;
    }

    //--------------------------------------------------------------------------
    // check for quick return
    //--------------------------------------------------------------------------

    if (nnz == 0)
    {
        (*component) = parent ;
        LAGraph_FREE_WORK ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // final phase
    //--------------------------------------------------------------------------

    bool changing = true ;
    while (changing)
    {
        // hooking & shortcutting
        // mngp = min (mngp, A*gp) using the MIN_SECOND semiring
        GrB_TRY (GrB_mxv (mngp, NULL, GrB_MIN_UINT64,
                          GrB_MIN_SECOND_SEMIRING_UINT64, A, gp, NULL)) ;

        // parent = min (parent, C*mngp) where C is C(i,j) = true if i=Px(j)
        GrB_TRY (Reduce_assign (parent, mngp, C, &Cp, &Px, &Cx, msg)) ;

        // parent = min (parent, mngp, gp)
        GrB_TRY (GrB_eWiseAdd (parent, NULL, GrB_MIN_UINT64, GrB_MIN_UINT64,
                               mngp, gp, NULL)) ;

        // calculate grandparent: gp_new = parent (parent)
        GrB_TRY (GrB_Vector_extractTuples (NULL, Px, &n, parent)) ;
        GrB_TRY (GrB_extract (gp_new, NULL, NULL, parent, Px, n, NULL)) ;

        // terminate if gp and gp_new are the same
        GrB_TRY (GrB_eWiseMult (mod, NULL, NULL, GrB_NE_UINT64, gp_new, gp,
            NULL)) ;
        GrB_TRY (GrB_reduce (&changing, NULL, GrB_LOR_MONOID_BOOL, mod, NULL)) ;

        // swap gp and gp_new
        GrB_Vector t = gp ; gp = gp_new ; gp_new = t ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*component) = parent ;
    LAGraph_FREE_WORK ;
    return (GrB_SUCCESS) ;
#endif
}
