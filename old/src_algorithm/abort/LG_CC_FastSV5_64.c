//------------------------------------------------------------------------------
// LG_CC_FastSV5_64: connected components (64-bit version)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Code is based on the algorithm described in the following paper
// Zhang, Azad, Hu. FastSV: FastSV: A Distributed-Memory Connected Component
// Algorithm with Fast Convergence (SIAM PP20)

// A subsequent update to the algorithm is here (which might not be reflected
// in this code):
//
// Yongzhe Zhang, Ariful Azad, Aydin Buluc: Parallel algorithms for finding
// connected components using linear algebra. J. Parallel Distributed Comput.
// 144: 14-27 (2020).

// Modified by Tim Davis, Texas A&M University

// The input graph G must be undirected, or directed but with an adjacency
// matrix with symmetric structure. Self-edges (diagonal entries) are OK, and
// are ignored.  The values and type of G->A are ignored; just its structure is
// accessed.

// This function cannot be called by multiple user threads, since it unpacks
// G->A and then packs it back.  G->A is unchanged when the function returns,
// but during execution G->A is empty.

#define LAGraph_FREE_ALL ;

#include "LG_internal.h"

#if !LG_VANILLA

#if (! LG_SUITESPARSE )
#error "SuiteSparse:GraphBLAS v6.0.0 or later required"
#endif

//------------------------------------------------------------------------------
// hash functions: todo describe me
//------------------------------------------------------------------------------

// hash table size must be a power of 2
#define HASH_SIZE 1024

// number of samples to insert into the hash table
// todo: this seems to be a lot of entries for a HASH_SIZE of 1024.
// There could be lots of collisions.
#define HASH_SAMPLES 864

#define HASH(x) (((x << 4) + x) & (HASH_SIZE-1))
#define NEXT(x) ((x + 23) & (HASH_SIZE-1))

//------------------------------------------------------------------------------
// ht_init: todo describe me
//------------------------------------------------------------------------------

// Clear the hash table counts (ht_val [0:HASH_SIZE-1] = 0), and set all hash
// table entries as empty (ht_key [0:HASH_SIZE-1] =-1).

// todo: the memset of ht_key is confusing

// todo: the name "ht_val" is confusing.  It is not a value, but a count of
// the number of times the value x = ht_key [h] has been inserted into the
// hth position in the hash table.  It should be renamed ht_cnt.

static inline void ht_init
(
    int64_t *ht_key,
    int64_t *ht_val
)
{
    memset (ht_key, -1, sizeof (int64_t) * HASH_SIZE) ;
    memset (ht_val,  0, sizeof (int64_t) * HASH_SIZE) ;
}

//------------------------------------------------------------------------------
// ht_sample: todo describe me
//------------------------------------------------------------------------------

//

static inline void ht_sample
(
    uint64_t *V,      // array of size n (todo: this is a bad variable name)
    int64_t n,
    int64_t samples,    // number of samples to take from V
    int64_t *ht_key,
    int64_t *ht_val,
    uint64_t *seed
)
{
    for (int64_t k = 0 ; k < samples ; k++)
    {
        // select an entry from V at random
        int64_t x = V [LAGraph_Random60 (seed) % n] ;

        // find x in the hash table
        int64_t h = HASH (x) ;
        while (ht_key [h] != -1 && ht_key [h] != x)
        {
            h = NEXT (h) ;
        }

        ht_key [h] = x ;
        ht_val [h]++ ;
    }
}

//------------------------------------------------------------------------------
// ht_most_frequent: todo describe me
//------------------------------------------------------------------------------

// todo what if key is returned as -1?  Code breaks.  todo: handle this case

static inline int64_t ht_most_frequent
(
    int64_t *ht_key,
    int64_t *ht_val
)
{
    int64_t key = -1 ;
    int64_t val = 0 ;                       // max (ht_val [0:HASH_SIZE-1])
    for (int64_t h = 0 ; h < HASH_SIZE ; h++)
    {
        if (ht_val [h] > val)
        {
            key = ht_key [h] ;
            val = ht_val [h] ;
        }
    }
    return (key) ;      // return most frequent key
}

//------------------------------------------------------------------------------
// Reduce_assign:  w (index) += s, using MIN as the "+=" accum operator
//------------------------------------------------------------------------------

// The index array, of size n can have duplicates.  The vectors w and s are
// full (all entries present).  This function computes:
//
//      for (j = 0 ; j < n ; j++)
//      {
//          uint64_t i = index [j] ;
//          w [i] = min (w [i], s [j]) ;
//      }
//
//  If C(i,j) = true where i == index [j], then this can be written with the
//  min_second semiring:
//
//      w = min (w, C*s)

static inline int Reduce_assign
(
    GrB_Vector w,           // vector of size n, all entries present
    GrB_Vector s,           // vector of size n, all entries present
    GrB_Matrix C,           // boolean matrix of size n-by-n
    GrB_Index **Cp_handle,  // array of size n+1, equal to 0:n
    GrB_Index **Ci_handle,  // index array of size n, can have duplicates
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
// LG_CC_FastSV5_64
//------------------------------------------------------------------------------

// The output of LG_CC_FastSV5 is a vector component, where
// component(i)=s if node i is in the connected compononent whose
// representative node is node s.  If s is a representative, then
// component(s)=s.  The number of connected components in the graph G is the
// number of representatives.

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL                            \
{                                                   \
    LAGraph_Free ((void **) &Cp) ;          \
    LAGraph_Free ((void **) &Cx) ;          \
    LAGraph_Free ((void **) &V) ;           \
    LAGraph_Free ((void **) &ht_key) ;      \
    LAGraph_Free ((void **) &ht_val) ;      \
    /* todo why is T not freed?? */                 \
    GrB_free (&t) ;                                 \
    GrB_free (&f) ;                                 \
    GrB_free (&gp) ;                                \
    GrB_free (&mngp) ;                              \
    GrB_free (&gp_new) ;                            \
    GrB_free (&mod) ;                               \
}

#endif

int LG_CC_FastSV5_64        // SuiteSparse:GraphBLAS method, with GxB extensions
(
    // output
    GrB_Vector *component,  // component(i)=s if node is in the component s
    // inputs
    LAGraph_Graph G,        // input graph, G->A can change
    char *msg
)
{

#if LG_VANILLA
    LG_CHECK (0, -1, "SuiteSparse required for this method") ;
#else

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    uint64_t *V = NULL ;
    int64_t *ht_key = NULL, *ht_val = NULL ;
    GrB_Index n, nnz ;
    GrB_Vector f = NULL, gp_new = NULL, mngp = NULL, mod = NULL, gp = NULL,
        t = NULL ;
    GrB_Matrix T = NULL, C = NULL ;
    GrB_Index *Cp = NULL ;
    GrB_Index Cp_size = 0 ;
    bool *Cx = NULL ;

    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;
    LG_CHECK (component == NULL, -1, "component parameter is NULL") ;

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
        LG_CHECK (false, -1, "input must be symmetric") ;
    }

    GrB_Matrix S = G->A ;
    GrB_TRY (GrB_Matrix_nrows (&n, S)) ;
    GrB_TRY (GrB_Matrix_nvals (&nnz, S)) ;

    #define FASTSV_SAMPLES 4

    bool sampling = (n * FASTSV_SAMPLES * 2 < nnz) ;

    // random number seed
    uint64_t seed = n ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // determine # of threads to use
    int nthreads ;
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads, NULL)) ;
    nthreads = LAGraph_MIN (nthreads, n / 16) ;
    nthreads = LAGraph_MAX (nthreads, 1) ;

    // vectors
    GrB_TRY (GrB_Vector_new (&f,      GrB_UINT64, n)) ;
    GrB_TRY (GrB_Vector_new (&gp_new, GrB_UINT64, n)) ;
    GrB_TRY (GrB_Vector_new (&mod,    GrB_BOOL,   n)) ;

    V = LAGraph_Malloc (n, sizeof (uint64_t)) ;
    GrB_TRY (GrB_assign (f, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_apply (f, NULL, NULL, GrB_ROWINDEX_INT64, f, 0, NULL)) ;
    GrB_TRY (GrB_Vector_extractTuples (NULL, V, &n, f)) ;
    bool V_is_identity = true ; // true if V is 0:n-1

    GrB_TRY (GrB_Vector_dup (&gp,   f)) ;
    GrB_TRY (GrB_Vector_dup (&mngp, f)) ;

    // allocate the hash table
    ht_key = LAGraph_Malloc (HASH_SIZE, sizeof (int64_t)) ;
    ht_val = LAGraph_Malloc (HASH_SIZE, sizeof (int64_t)) ;
    LG_CHECK (ht_key == NULL || ht_val == NULL, -1, "out of memory") ;

    // create Cp = 0:n, and Cx = true, and the empty C matrix
    GrB_TRY (GrB_Vector_new (&t, GrB_INT64, n+1)) ;
    GrB_TRY (GrB_assign (t, NULL, NULL, 0, GrB_ALL, n+1, NULL)) ;
    GrB_TRY (GrB_apply (t, NULL, NULL, GrB_ROWINDEX_INT64, t, 0, NULL)) ;
    GrB_TRY (GxB_Vector_unpack_Full (t, (void **) &Cp, &Cp_size, NULL, NULL)) ;
    Cx = (bool *) LAGraph_Malloc (1, sizeof (bool)) ;
    Cx [0] = true ;
    GrB_TRY (GrB_free (&t)) ;
    GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;

    //--------------------------------------------------------------------------
    // sample phase
    //--------------------------------------------------------------------------

    if (sampling)
    {

        //----------------------------------------------------------------------
        // export S = G->A in CSR format
        //----------------------------------------------------------------------

        // S is not modified.  It is only exported so that its contents can be
        // read by the parallel loops below.

        GrB_Type type ;
        GrB_Index nrows, ncols, nvals ;
        size_t typesize ;
        int64_t nonempty ;
        GrB_Index *Sp, *Sj ;
        void *Sx ;
        bool S_jumbled = false ;
        GrB_Index Sp_size, Sj_size, Sx_size ;
        bool S_iso = false ;

        GrB_TRY (GrB_Matrix_nvals (&nvals, S)) ;
        GrB_TRY (GxB_Matrix_export_CSR (&S, &type, &nrows, &ncols, &Sp, &Sj,
            &Sx, &Sp_size, &Sj_size, &Sx_size,
            &S_iso, &S_jumbled, NULL)) ;
        GrB_TRY (GxB_Type_size (&typesize, type)) ;
        G->A = NULL ;

        //----------------------------------------------------------------------
        // allocate space to construct T
        //----------------------------------------------------------------------

        GrB_Index Tp_len = nrows+1, Tp_size = Tp_len*sizeof(GrB_Index);
        GrB_Index Tj_len = nvals,   Tj_size = Tj_len*sizeof(GrB_Index);
        GrB_Index Tx_len = nvals ;
        GrB_Index *Tp = LAGraph_Malloc (Tp_len, sizeof (GrB_Index)) ;
        GrB_Index *Tj = LAGraph_Malloc (Tj_len, sizeof (GrB_Index)) ;
        GrB_Index Tx_size = typesize ;
        void *Tx = LAGraph_Calloc (1, typesize) ;   // T is iso

        // todo check out-of-memory conditions

        //----------------------------------------------------------------------
        // allocate workspace
        //----------------------------------------------------------------------

        int64_t *range = LAGraph_Malloc (nthreads + 1, sizeof (int64_t)) ;
        GrB_Index *count = LAGraph_Malloc (nthreads + 1, sizeof (GrB_Index)) ;
        // todo check out-of-memory conditions

        memset (count, 0, sizeof (GrB_Index) * (nthreads + 1)) ;

        //----------------------------------------------------------------------
        // define parallel tasks to construct T
        //----------------------------------------------------------------------

        // thread tid works on rows range[tid]:range[tid+1]-1 of S and T
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

        // T (i,:) consists of the first FASTSV_SAMPLES of S (i,:).

        // todo: this could be done by GxB_Select, using a new operator.  Need
        // to define a set of GxB_SelectOp operators that would allow for this.

        // Note that Tx is not modified.  Only Tp and Tj are constructed.

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = count [tid] ;
            Tp [range [tid]] = p ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                // construct T (i,:) from the first entries in S (i,:)
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

        // Note that Tx is unmodified.

        // in SuiteSparse:GraphBLAS v5, sizes are in bytes, not entries
        GrB_Index Tp_siz = Tp_size ;
        GrB_Index Tj_siz = Tj_size ;
        GrB_Index Tx_siz = Tx_size ;

        GrB_Index t_nvals = Tp [nrows] ;
        GrB_TRY (GxB_Matrix_import_CSR (&T, type, nrows, ncols,
                &Tp, &Tj, &Tx, Tp_siz, Tj_siz, Tx_siz,
                true,   // T is iso
                S_jumbled, NULL)) ;

        //----------------------------------------------------------------------
        // find the connected components of T
        //----------------------------------------------------------------------

        // TODO: this is identical to the final phase below. Make it a function

        bool diff = true ;
        while (diff)
        {
            // hooking & shortcutting
            // mngp = min (mngp, T*gp) using the MIN_SECOND semiring
            GrB_TRY (GrB_mxv (mngp, NULL, GrB_MIN_UINT64,
                GrB_MIN_SECOND_SEMIRING_UINT64, T, gp, NULL)) ;
            if (!V_is_identity)
            {
                // f = min (f, C*mngp) where C is C(i,j) = true if i=V(j)
                LAGraph_TRY (Reduce_assign (f, mngp, C, &Cp, &V, &Cx, msg)) ;
            }

            // f = min (f, mngp, gp)
            GrB_TRY (GrB_eWiseAdd (f, NULL, GrB_MIN_UINT64, GrB_MIN_UINT64,
                mngp, gp, NULL)) ;

            // calculate grandparent: gp_new = f (f)
            GrB_TRY (GrB_Vector_extractTuples (NULL, V, &n, f)) ;
            GrB_TRY (GrB_extract (gp_new, NULL, NULL, f, V, n, NULL)) ;
            V_is_identity = false ;

            // terminate if gp and gb_new are the same
            GrB_TRY (GrB_eWiseMult (mod, NULL, NULL, GrB_NE_UINT64, gp_new,
                gp, NULL)) ;
            GrB_TRY (GrB_reduce (&diff, NULL, GrB_LOR_MONOID_BOOL, mod,
                NULL)) ;

            // swap gp and gp_new
            GrB_Vector t = gp ; gp = gp_new ; gp_new = t ;
        }

        //----------------------------------------------------------------------
        // estimate the largest connected component
        //----------------------------------------------------------------------

        // key is set to the representative node of the largest connected
        // component.  Since sampling is used, this is an estimate

        ht_init (ht_key, ht_val) ;
        ht_sample (V, n, HASH_SAMPLES, ht_key, ht_val, &seed) ;
        int64_t key = ht_most_frequent (ht_key, ht_val) ;
        // todo: what if key is returned as -1?  Then T below is invalid.

        //----------------------------------------------------------------------
        // collapse the largest connected component
        //----------------------------------------------------------------------

        // All edges in, or to, the largest connected component is removed from
        // T.  Next, if the row T(i,:) has enough space, and node i is not in
        // the largest connected component, a single edge T(i,key) = true is
        // added.

        int64_t t_nonempty = -1 ;
        bool T_jumbled = false, T_iso = true ;

        // export T
        GrB_TRY (GxB_Matrix_export_CSR (&T, &type, &nrows, &ncols, &Tp, &Tj,
            &Tx, &Tp_siz, &Tj_siz, &Tx_siz,
            &T_iso, &T_jumbled, NULL)) ;

        // FIXME: This parallel loop is badly load balanced.  Each thread
        // operates on the same number of rows of S, regardless of how many
        // entries appear in each set of rows.  It uses one thread per task,
        // statically scheduled.

        #pragma omp parallel for num_threads(nthreads) schedule(static) \
            reduction(||:T_jumbled)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index ptr = Sp [range [tid]] ;
            // thread tid scans S (range [tid]:range [tid+1]-1,:),
            // and constructs T(i,:) for all rows in this range.
            bool my_T_jumbled = false ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                int64_t pv = V [i] ;  // what is pv?
                Tp [i] = ptr ;          // start the construction of T(i,:)
                // T(i,:) is empty if pv == key
                if (pv != key)
                {
                    // scan S(i,:)
                    for (GrB_Index p = Sp [i] ; p < Sp [i+1] ; p++)
                    {
                        // get S(i,j)
                        int64_t j = Sj [p] ;
                        if (V [j] != key)
                        {
                            // add the entry T(i,j) to T, but skip it if
                            // V [j] is equal to key
                            Tj [ptr++] = j ;
                        }
                    }
                    // add the entry T(i,key) if there is room for it in T(i,:)
                    if (ptr - Tp [i] < Sp [i+1] - Sp [i])
                    {
                        Tj [ptr++] = key ;
                        // this step can cause T to become jumbled
                        my_T_jumbled = true ;
                    }
                }
            }
            T_jumbled = T_jumbled || my_T_jumbled ;
            // count the number of entries inserted into T by this thread
            count [tid] = ptr - Tp [range [tid]] ;
        }

        // Compact empty space out of Tj not filled in from the above phase.
        // This is a lot of work and should be done in parallel.
        GrB_Index offset = 0 ;
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            memcpy (Tj + offset, Tj + Tp [range [tid]],
                sizeof (GrB_Index) * count [tid]) ;
            offset += count [tid] ;
            count [tid] = offset - count [tid] ;
        }

        // Compact empty space out of Tp
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index ptr = Tp [range [tid]] ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                Tp [i] -= ptr - count [tid] ;
            }
        }

        // finalize T
        Tp [n] = offset ;

        // free workspace
        LAGraph_Free ((void **) &count) ;
        LAGraph_Free ((void **) &range) ;

        // import S (unchanged since last export)
        GrB_TRY (GxB_Matrix_import_CSR (&S, type, nrows, ncols,
                &Sp, &Sj, &Sx, Sp_size, Sj_size, Sx_size,
                S_iso, S_jumbled, NULL)) ;

        // import T for the final phase
        GrB_TRY (GxB_Matrix_import_CSR (&T, type, nrows, ncols,
                &Tp, &Tj, &Tx, Tp_siz, Tj_siz, Tx_siz,
                T_iso, T_jumbled, NULL)) ;

        // restore G->A
        G->A = S ;

    }
    else
    {

        // no sampling; the final phase operates on the whole graph
        T = S ;

    }

    //--------------------------------------------------------------------------
    // final phase
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_nvals (&nnz, T)) ;

    bool diff = true ;
    while (diff && nnz > 0)
    {
        // hooking & shortcutting
        // mngp = min (mngp, T*gp) using the MIN_SECOND semiring
        GrB_TRY (GrB_mxv (mngp, NULL, GrB_MIN_UINT64,
                          GrB_MIN_SECOND_SEMIRING_UINT64, T, gp, NULL)) ;

        if (!V_is_identity)
        {
            // f = min (f, C*mngp) where C is C(i,j) = true if i=V(j)
            GrB_TRY (Reduce_assign (f, mngp, C, &Cp, &V, &Cx, msg)) ;
            V_is_identity = false ;
        }

        // f = min (f, mngp, gp)
        GrB_TRY (GrB_eWiseAdd (f, NULL, GrB_MIN_UINT64, GrB_MIN_UINT64,
                               mngp, gp, NULL)) ;

        // calculate grandparent: gp_new = f (f)
        GrB_TRY (GrB_Vector_extractTuples (NULL, V, &n, f)) ;
        GrB_TRY (GrB_extract (gp_new, NULL, NULL, f, V, n, NULL)) ;

        // terminate if gp and gb_new are the same
        GrB_TRY (GrB_eWiseMult (mod, NULL, NULL, GrB_NE_UINT64, gp_new, gp,
            NULL)) ;
        GrB_TRY (GrB_reduce (&diff, NULL, GrB_LOR_MONOID_BOOL, mod, NULL)) ;

        // swap gp and gp_new
        GrB_Vector t = gp ; gp = gp_new ; gp_new = t ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*component) = f ;
    f = NULL ;
    if (sampling)
    {
        GrB_free (&T) ;
    }
    LAGraph_FREE_ALL ;
    return (0) ;
#endif
}
