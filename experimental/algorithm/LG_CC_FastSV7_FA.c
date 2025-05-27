//------------------------------------------------------------------------------
// LG_CC_FastSV7: connected components
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Yongzhe Zhang, modified by Timothy A. Davis, Texas A&M
// University

//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->is_symmetric_structure must be known),
// but it is not user-callable (see LAGr_ConnectedComponents instead).

// Code is based on the algorithm described in the following paper:
// Zhang, Azad, Hu. FastSV: A Distributed-Memory Connected Component
// Algorithm with Fast Convergence (SIAM PP20)

// A subsequent update to the algorithm is here (which might not be reflected
// in this code):
// Yongzhe Zhang, Ariful Azad, Aydin Buluc: Parallel algorithms for finding
// connected components using linear algebra. J. Parallel Distributed Comput.
// 144: 14-27 (2020).

// Modified by Tim Davis, Texas A&M University: revised Reduce_assign to use
// purely GrB* and GxB* methods and the matrix Parent.  Added warmup phase.
// Changed to use GxB load/unload.  Converted to use the LAGraph_Graph object.
// Exploiting iso status for the temporary matrices Parent and T.

// Modified by Gabriel Gomez, Texas A&M University: moved Parent matrix trick 
// out to LAGraph_FastAssign.

// The input graph G must be undirected, or directed and with an adjacency
// matrix that has a symmetric structure.  Self-edges (diagonal entries) are
// OK, and are ignored.  The values and type of A are ignored; just its
// structure is accessed.

// NOTE: This function must not be called by multiple user threads at the same
// time on the same graph G, since it unloads G->A and loads it back when
// done.  G->A is unchanged when the function returns, but during execution
// G->A is empty.  This will be fixed once the todos are finished below, and
// G->A will then become a truly read-only object (assuming GrB_wait (G->A)
// has been done first).

// #define TIMINGS

#define __STDC_WANT_LIB_EXT1__ 1
#include <string.h>

#define LG_FREE_ALL ;
#include "LG_internal.h"
#include "LAGraphX.h"

static double timings [16] ;

#if LG_SUITESPARSE_GRAPHBLAS_V10

//==============================================================================
// fastsv: find the components of a graph
//==============================================================================

static inline GrB_Info fastsv
(
    GrB_Matrix A,           // adjacency matrix, G->A or a subset of G->A
    GrB_Vector parent2,     // workspace
    GrB_Vector mngp,        // min neighbor grandparent
    GrB_Vector *gp,         // grandparent
    GrB_Vector *gp_new,     // new grandparent (swapped with gp)
    GrB_Vector t,           // workspace
    GrB_BinaryOp eq,        // GrB_EQ_(integer type)
    GrB_BinaryOp min,       // GrB_MIN_(integer type)
    GrB_Semiring min_2nd,   // GrB_MIN_SECOND_(integer type)
    GrB_Vector parent,      // parent
    GrB_Vector ramp,        // [0:n] used to speed up FastAssign
    char *msg
)
{
    bool done = false ;
//  #ifdef TIMINGS
//  int pass = 0 ;
//  #endif

    while (true)
    {
//      #ifdef TIMINGS
//      printf ("\n-------------------------------------------fastsv: %d\n",
//          ++pass) ;
//      #endif

        //----------------------------------------------------------------------
        // hooking & shortcutting
        //----------------------------------------------------------------------

        // mngp = min (mngp, A*gp) using the MIN_SECOND semiring
        GRB_TRY (GrB_mxv (mngp, NULL, min, min_2nd, A, *gp, NULL)) ;

        //----------------------------------------------------------------------
        // parent2 = min (mngp, gp)
        //----------------------------------------------------------------------

        // The parent vector should not be allised into FastAssign, so the 
        // accumulation is done in a workspace vector, parent2.

        GRB_TRY (GrB_eWiseAdd (parent2, NULL, NULL, min, mngp, *gp, NULL)) ;

        // LAGraph_FastAssign: This function computes the following, which
        // is done explicitly in the Reduce_assign function in LG_CC_Boruvka:
        //
        //      for (j = 0 ; j < n ; j++)
        //      {
        //          uint64_t i = parent [j] ;
        //          parent2 [i] = min (parent2 [i], mngp [j]) ;
        //      }
        //
        // LAGraph_FastAssign does this by building a matrix. 
        // (See LAGraph_FastAssign.c) 
        // Giving it a full ramp vector speeds up the function

        LG_TRY (LAGraph_FastAssign_Semiring(
            parent2, NULL, min, parent, mngp, ramp, min_2nd, NULL, msg));

        //----------------------------------------------------------------------
        // parent = min (parent, parent2)
        //----------------------------------------------------------------------

        GRB_TRY (GrB_assign (parent, NULL, min, parent2, GrB_ALL, 0, NULL)) ;

        //----------------------------------------------------------------------
        // calculate grandparent: gp_new = parent (parent)
        //----------------------------------------------------------------------

        GRB_TRY (GrB_extract (*gp_new, NULL, NULL, parent, parent, NULL)) ;

        //----------------------------------------------------------------------
        // terminate if gp and gp_new are the same
        //----------------------------------------------------------------------
        GRB_TRY (GrB_eWiseMult (t, NULL, NULL, eq, *gp_new, *gp, NULL)) ;
        GRB_TRY (GrB_reduce (&done, NULL, GrB_LAND_MONOID_BOOL, t, NULL)) ;
        if (done) break ;

        // swap gp and gp_new
        GrB_Vector s = (*gp) ; (*gp) = (*gp_new) ; (*gp_new) = s ;
    }
    return (GrB_SUCCESS) ;
}

//==============================================================================
// LG_CC_FastSV7
//==============================================================================

// The output of LG_CC_FastSV* is a vector component, where component(i)=r if
// node i is in the connected compononent whose representative is node r.  If r
// is a representative, then component(r)=r.  The number of connected
// components in the graph G is the number of representatives.

#undef  LG_FREE_WORK
#define LG_FREE_WORK                            \
{                                               \
    LAGraph_Free ((void **) &Tp, NULL) ;        \
    LAGraph_Free ((void **) &Tj, NULL) ;        \
    LAGraph_Free ((void **) &Tx, NULL) ;        \
    LAGraph_Free ((void **) &ht_key, NULL) ;    \
    LAGraph_Free ((void **) &ht_count, NULL) ;  \
    LAGraph_Free ((void **) &count, NULL) ;     \
    LAGraph_Free ((void **) &range, NULL) ;     \
    LAGraph_Free ((void **) &Px, NULL) ;        \
    GrB_free (&T) ;                             \
    GrB_free (&t) ;                             \
    GrB_free (&gp) ;                            \
    GrB_free (&mngp) ;                          \
    GrB_free (&gp_new) ;                        \
    GrB_free (&parent2) ;                       \
    GrB_free (&ramp_v) ;                        \
    GrB_free (&A_Container) ;                   \
    GrB_free (&T_Container) ;                   \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&parent) ;                        \
    LG_FREE_WORK ;                              \
}

#endif

// get/set macros for 32/64 bit arrays:
#define AP(k) (Ap32 ? Ap32 [k] : Ap64 [k])
#define AJ(k) (Aj32 ? Aj32 [k] : Aj64 [k])
#define PARENT(i) (Px32 ? Px32 [i] : Px64 [i])
#define TP(k) (Tp32 ? Tp32 [k] : Tp64 [k])
#define TJ(k) (Tj32 ? Tj32 [k] : Tj64 [k])
#define SET_TP(k,p) { if (Tp32) { Tp32 [k] = p ; } else { Tp64 [k] = p ; }}
#define SET_TJ(k,i) { if (Tj32) { Tj32 [k] = i ; } else { Tj64 [k] = i ; }}

#ifdef TIMINGS
static void print_timings (double timings [16])
{
    double total = timings [0] + timings [1] + timings [2] ;
    printf ("SV7 %12.6f (%4.1f%%) init\n", timings [0], 100. * timings [0] / total) ;
    printf ("SV7 %12.6f (%4.1f%%) total sampling:\n", timings [1], 100. * timings [1] / total) ;
    printf ("SV7        %12.6f (%4.1f%%) setup T\n", timings [3], 100. * timings [3] / total) ;
    printf ("SV7        %12.6f (%4.1f%%) create T\n", timings [4], 100. * timings [4] / total) ;
    printf ("SV7        %12.6f (%4.1f%%) fastsv sample\n", timings [5], 100 * timings [5] / total) ;
    printf ("SV7        %12.6f (%4.1f%%) hash\n", timings [6], 100. * timings [6] / total) ;
    printf ("SV7        %12.6f (%4.1f%%) prune\n", timings [7], 100. * timings [7] / total) ;
    printf ("SV7 %12.6f (%4.1f%%) total final\n", timings [2], 100. * timings [2] / total) ;
}
#endif

int LG_CC_FastSV7_FA         // SuiteSparse:GraphBLAS method, with GraphBLAS v10
(
    // output:
    GrB_Vector *component,  // component(i)=r if node is in the component r
    // input:
    LAGraph_Graph G,        // input graph (modified then restored)
    char *msg
)
{

#if LG_SUITESPARSE_GRAPHBLAS_V10

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    #ifdef TIMINGS
    double timings [16] ;
    for (int kk = 0 ; kk < 16 ; kk++) timings [kk] = 0 ;
    double tic = LAGraph_WallClockTime ( ) ;
    LG_SET_BURBLE (false) ;
    #endif

    int64_t *range = NULL ;
    void *Px = NULL ;
    uint64_t Px_size = 0 ;
    GrB_Index n, nvals, *ht_key = NULL, *count = NULL ;
    void *Tp = NULL, *Tj = NULL ;
    GrB_Vector parent = NULL, gp_new = NULL, mngp = NULL, gp = NULL, t = NULL,
        parent2 = NULL, ramp_v = NULL ;
    GrB_Matrix T = NULL ;
    void *Tx = NULL ;
    int *ht_count = NULL ;
    GxB_Container A_Container = NULL, T_Container = NULL ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (component != NULL, GrB_NULL_POINTER) ;
    (*component) = NULL ;

    LG_ASSERT_MSG ((G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->is_symmetric_structure == LAGraph_TRUE)),
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be known to be symmetric") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;

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
        // use 64-bit integers
        Uint = GrB_UINT64 ;
        Int  = GrB_INT64  ;
        ramp = GrB_ROWINDEX_INT64 ;
        min  = GrB_MIN_INT64 ;
        imin = GrB_MIN_INT64 ;
        eq   = GrB_EQ_INT64 ;
        min_2nd  = GrB_MIN_SECOND_SEMIRING_INT64 ;
        min_2ndi = GxB_MIN_SECONDI_INT64 ;
    }
    else
    {
        // use 32-bit integers
        Uint = GrB_UINT32 ;
        Int  = GrB_INT32  ;
        ramp = GrB_ROWINDEX_INT32 ;
        min  = GrB_MIN_INT32 ;
        imin = GrB_MIN_INT32 ;
        eq   = GrB_EQ_INT32 ;
        min_2nd  = GrB_MIN_SECOND_SEMIRING_INT32 ;
        min_2ndi = GxB_MIN_SECONDI_INT32 ;
    }

    // FASTSV_SAMPLES: number of samples to take from each row A(i,:).
    // Sampling is used if the average degree is > 8 and if n > 1024.
    #define FASTSV_SAMPLES 4
    bool sampling = (nvals > n * FASTSV_SAMPLES * 2 && n > 1024) ;

    //--------------------------------------------------------------------------
    // make ramp needed for FastAssign speedup
    //--------------------------------------------------------------------------
    GRB_TRY (GrB_Vector_new (&(ramp_v), Uint, n+1)) ;
    GRB_TRY (GrB_assign (ramp_v, NULL, NULL, 0, GrB_ALL, n+1,
        NULL)) ;
    GRB_TRY (GrB_apply (ramp_v, NULL, NULL, ramp, ramp_v, 0, NULL)) ;
// [ todo: nthreads will not be needed once GxB_select with a GxB_RankUnaryOp
// and a new GxB_extract are added to SuiteSparse:GraphBLAS.
    // determine # of threads to use
    int nthreads, nthreads_outer, nthreads_inner ;
    LG_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
    nthreads = nthreads_outer * nthreads_inner ;
    nthreads = LAGRAPH_MIN (nthreads, n / 16) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;
// ]

    GRB_TRY (GxB_Container_new (&A_Container)) ;
    GRB_TRY (GxB_Container_new (&T_Container)) ;

    //--------------------------------------------------------------------------
    // warmup: parent = min (0:n-1, A*1) using the MIN_SECONDI semiring
    //--------------------------------------------------------------------------

    // parent (i) = min (i, j) for all entries A(i,j).  This warmup phase takes only
    // O(n) time, because of how the MIN_SECONDI semiring is implemented in
    // SuiteSparse:GraphBLAS.  A is held by row, and the first entry in A(i,:)
    // is the minimum index j, so only the first entry in A(i,:) needs to be
    // considered for each row i.

    GRB_TRY (GrB_Vector_new (&t, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&parent, Int, n)) ;
    GRB_TRY (GrB_assign (t, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_assign (parent, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_apply (parent, NULL, NULL, ramp, parent, 0, NULL)) ;
    GRB_TRY (GrB_mxv (parent, NULL, imin, min_2ndi, A, t, NULL)) ;
    GRB_TRY (GrB_free (&t)) ;

    // copy parent into gp and mngp.
    GRB_TRY (GrB_Vector_dup (&gp, parent)) ;
    GRB_TRY (GrB_Vector_dup (&mngp, parent)) ;

    // allocate workspace vectors
    GRB_TRY (GrB_Vector_new (&gp_new, Int, n)) ;
    GRB_TRY (GrB_Vector_new (&t, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&parent2, Int, n)) ;

    #ifdef TIMINGS
    double toc = LAGraph_WallClockTime ( ) ;
    timings [0] = toc - tic ;  // init time
    tic = toc ;
    #endif

    //--------------------------------------------------------------------------
    // sample phase
    //--------------------------------------------------------------------------

    if (sampling)
    {

// [ todo: GxB_select, using a new operator: GxB_RankUnaryOp, will do all this,
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

        //----------------------------------------------------------------------
        // unload A in CSR format
        //----------------------------------------------------------------------

        #ifdef TIMINGS
        double tic2 = LAGraph_WallClockTime ( ) ;
        #endif

        void *Ap = NULL, *Aj = NULL ;
        uint64_t Ap_size, Aj_size, Ap_len, Aj_len ;
        bool A_jumbled, A_iso ;
        int Ap_handling, Aj_handling ;

        // unload A in sparse CSR format into the A_Container
        GRB_TRY (GrB_set (A, GxB_SPARSE, GxB_SPARSITY_CONTROL)) ;
        GRB_TRY (GrB_set (A, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
        GRB_TRY (GxB_unload_Matrix_into_Container (A, A_Container, NULL)) ;
        A_jumbled = A_Container->jumbled ;
        A_iso = A_Container->iso ;

        // unload A_Container->p,i into the C arrays, Ap and Aj
        GrB_Type Ap_type, Aj_type ;
        GRB_TRY (GxB_Vector_unload (A_Container->p, &Ap, &Ap_type, &Ap_len,
            &Ap_size, &Ap_handling, NULL)) ;
        GRB_TRY (GxB_Vector_unload (A_Container->i, &Aj, &Aj_type, &Aj_len,
            &Aj_size, &Aj_handling, NULL)) ;

        bool Ap_is_32 = (Ap_type == GrB_UINT32 || Ap_type == GrB_INT32) ;
        bool Aj_is_32 = (Aj_type == GrB_UINT32 || Aj_type == GrB_INT32) ;
        const uint32_t *Ap32 = Ap_is_32 ? Ap : NULL ;
        const uint64_t *Ap64 = Ap_is_32 ? NULL : Ap ;
        const uint32_t *Aj32 = Aj_is_32 ? Aj : NULL ;
        const uint64_t *Aj64 = Aj_is_32 ? NULL : Aj ;

        //----------------------------------------------------------------------
        // allocate workspace, including space to construct T
        //----------------------------------------------------------------------

        bool Tp_is_32 = (nvals < UINT32_MAX) ;
        bool Tj_is_32 = (n < INT32_MAX) ;
        GrB_Type Tp_type = Tp_is_32 ? GrB_UINT32 : GrB_UINT64 ;
        GrB_Type Tj_type = Tj_is_32 ? GrB_UINT32 : GrB_UINT64 ;
        size_t tpsize = Tp_is_32 ? sizeof (uint32_t) : sizeof (uint64_t) ;
        size_t tjsize = Tj_is_32 ? sizeof (uint32_t) : sizeof (uint64_t) ;

        GrB_Index Tp_size = (n+1) * tpsize ;
        GrB_Index Tj_size = nvals * tjsize ;
        GrB_Index Tx_size = sizeof (bool) ;
        LG_TRY (LAGraph_Malloc ((void **) &Tp, n+1, tpsize, msg)) ;
        LG_TRY (LAGraph_Malloc ((void **) &Tj, nvals, tjsize, msg)) ;
        LG_TRY (LAGraph_Calloc ((void **) &Tx, 1, sizeof (bool), msg)) ;
        LG_TRY (LAGraph_Malloc ((void **) &range, nthreads+1, sizeof (int64_t),
            msg)) ;
        LG_TRY (LAGraph_Calloc ((void **) &count, nthreads+1, sizeof (uint64_t),
            msg)) ;

        uint32_t *Tp32 = Tp_is_32 ? Tp : NULL ;
        uint64_t *Tp64 = Tp_is_32 ? NULL : Tp ;
        uint32_t *Tj32 = Tj_is_32 ? Tj : NULL ;
        uint64_t *Tj64 = Tj_is_32 ? NULL : Tj ;

        //----------------------------------------------------------------------
        // define parallel tasks to construct T
        //----------------------------------------------------------------------

        // thread tid works on rows range[tid]:range[tid+1]-1 of A and T
        int tid;
        for (tid = 0 ; tid <= nthreads ; tid++)
        {
            range [tid] = (n * tid + nthreads - 1) / nthreads ;
        }

        //----------------------------------------------------------------------
        // determine the number entries to be constructed in T for each thread
        //----------------------------------------------------------------------

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (tid = 0 ; tid < nthreads ; tid++)
        {
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                int64_t deg = AP (i + 1) - AP (i) ;
                count [tid + 1] += LAGRAPH_MIN (FASTSV_SAMPLES, deg) ;
            }
        }

        //----------------------------------------------------------------------
        // count = cumsum (count)
        //----------------------------------------------------------------------

        for (tid = 0 ; tid < nthreads ; tid++)
        {
            count [tid + 1] += count [tid] ;
        }

        #ifdef TIMINGS
        double toc2 = LAGraph_WallClockTime ( ) ;
        timings [3] = toc2 - tic2 ;  // setup T
        tic2 = toc2 ;
        #endif

        //----------------------------------------------------------------------
        // construct T
        //----------------------------------------------------------------------

        // T (i,:) consists of the first FASTSV_SAMPLES of A (i,:).

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (tid = 0 ; tid < nthreads ; tid++)
        {
            GrB_Index p = count [tid] ;
            int64_t ktid = range [tid] ;
            SET_TP (ktid, p) ;      // Tp [ktid] = p ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                // construct T (i,:) from the first entries in A (i,:)
                for (int64_t j = 0 ;
                    j < FASTSV_SAMPLES && AP (i) + j < AP (i + 1) ; j++)
                {
                    uint64_t pi = AP (i) + j ;
                    uint64_t j = AJ (pi) ;
                    SET_TJ (p, j) ;         // Tj [p] = j ;
                    p++ ;
                }
                SET_TP (i+1, p) ;           // Tp [i + 1] = p ;
            }
        }

        //----------------------------------------------------------------------
        // import the result into the GrB_Matrix T
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;

        uint64_t T_nvals = TP (n) ;
        uint64_t Tp_len = n+1 ;
        uint64_t Tj_len = T_nvals ;

        T_Container->nrows = n ;
        T_Container->ncols = n ;
        T_Container->nrows_nonempty = -1 ;
        T_Container->ncols_nonempty = -1 ;
        T_Container->nvals = T_nvals ;
        T_Container->format = GxB_SPARSE ;
        T_Container->orientation = GrB_ROWMAJOR ;
        T_Container->iso = true ;
        T_Container->jumbled = A_jumbled ;

        // load Tp, Tj, and Tx into the T_Container
        GRB_TRY (GxB_Vector_load (T_Container->p, &Tp, Tp_type, Tp_len,
            Tp_size, GrB_DEFAULT, NULL)) ;
        GRB_TRY (GxB_Vector_load (T_Container->i, &Tj, Tj_type, Tj_len,
            Tj_size, GrB_DEFAULT, NULL)) ;
        GRB_TRY (GxB_Vector_load (T_Container->x, &Tx, GrB_BOOL, 1,
            Tx_size, GrB_DEFAULT, NULL)) ;

        // load T from the T_Container
        GRB_TRY (GxB_load_Matrix_from_Container (T, T_Container, NULL)) ;

// ] todo: the above will all be done as a single call to GxB_select.

        #ifdef TIMINGS
        toc2 = LAGraph_WallClockTime ( ) ;
        timings [4] = toc2 - tic2 ;  // create T
        tic2 = toc2 ;
        #endif

        //----------------------------------------------------------------------
        // find the connected components of T
        //----------------------------------------------------------------------

        LG_TRY (fastsv (T, parent2, mngp, &gp, &gp_new, t, eq, min, min_2nd,
            parent, ramp_v, msg)) ;

        #ifdef TIMINGS
        toc2 = LAGraph_WallClockTime ( ) ;
        timings [5] = toc2 - tic2 ;  // fastsv, in sampling
        tic2 = toc2 ;
        #endif

        //----------------------------------------------------------------------
        // unload the parent i vector into the Px array
        //----------------------------------------------------------------------

        int handling = 0 ;
        GrB_Type type = NULL ;
        GRB_TRY (GxB_Vector_unload (parent, &Px, &type, &n, &Px_size,
            &handling, NULL)) ;
        bool Px_is_32 = (type == GrB_UINT32 || type == GrB_INT32) ;
        uint32_t *Px32 = Px_is_32 ? Px : NULL ;
        uint64_t *Px64 = Px_is_32 ? NULL : Px ;

        // At this point, the Px array holds the content of parent vector.

        //----------------------------------------------------------------------
        // use sampling to estimate the largest connected component in T
        //----------------------------------------------------------------------

        // The sampling below computes an estimate of the mode of the parent
        // vector, the contents of which are currently in the non-opaque Px
        // array.

        // hash table size must be a power of 2
        #define HASH_SIZE 1024
        // number of samples to insert into the hash table
        #define HASH_SAMPLES 864
        #define HASH(x) (((x << 4) + x) & (HASH_SIZE-1))
        #define NEXT(x) ((x + 23) & (HASH_SIZE-1))

        // allocate and initialize the hash table
        LG_TRY (LAGraph_Malloc ((void **) &ht_key, HASH_SIZE,
            sizeof (GrB_Index), msg)) ;
        LG_TRY (LAGraph_Calloc ((void **) &ht_count, HASH_SIZE,
            sizeof (int), msg)) ;
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
            // select an entry ii from PARENT at random
            uint64_t i = LG_Random64 (&seed) % n ;
            GrB_Index x = PARENT (i) ;
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

        #ifdef TIMINGS
        toc2 = LAGraph_WallClockTime ( ) ;
        timings [6] = toc2 - tic2 ;  // hash
        tic2 = toc2 ;
        #endif

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

// [ todo: replace this with GxB_extract with GrB_Vector index arrays.
// See https://github.com/GraphBLAS/graphblas-api-c/issues/67 .
// This method will not insert the new entries T(i,key) for rows i that have
// had entries deleted.  That can be done with GrB_assign, with an n-by-1 mask
// M computed from the before-and-after row degrees of A and T:
// M = (parent != key) && (out_degree(T) < out_degree(A))
// J [0] = key.
// GxB_Matrix_subassign_BOOL (T, M, NULL, true, GrB_ALL, n, J, 1, NULL)
// or with
// GrB_Col_assign (T, M, NULL, t, GrB_ALL, j, NULL) with an all-true
// vector t.

        // unload T from the T_Container; its contents are revised below
        GRB_TRY (GxB_unload_Matrix_into_Container (T, T_Container, NULL)) ;

        // unload Tp and Tj from the T_Container
        int ignore ;
        GRB_TRY (GxB_Vector_unload (T_Container->p, &Tp, &Tp_type, &Tp_len,
            &Tp_size, &ignore, NULL)) ;
        GRB_TRY (GxB_Vector_unload (T_Container->i, &Tj, &Tj_type, &Tj_len,
            &Tj_size, &ignore, NULL)) ;

        // these are likely to be unchanged since the last load of T
        Tp_is_32 = (Tp_type == GrB_UINT32 || Tp_type == GrB_INT32) ;
        Tj_is_32 = (Tj_type == GrB_UINT32 || Tj_type == GrB_INT32) ;
        Tp32 = Tp_is_32 ? Tp : NULL ;
        Tp64 = Tp_is_32 ? NULL : Tp ;
        Tj32 = Tj_is_32 ? Tj : NULL ;
        Tj64 = Tj_is_32 ? NULL : Tj ;

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (tid = 0 ; tid < nthreads ; tid++)
        {
            uint64_t ktid = range [tid] ;
            GrB_Index p = AP (ktid) ;
            // thread tid scans A (range [tid]:range [tid+1]-1,:),
            // and constructs T(i,:) for all rows in this range.
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                int64_t pi = PARENT (i) ;
                int64_t pstart = p ;
                SET_TP (i, p) ; // Tp [i] = p ; start the construction of T(i,:)
                // T(i,:) is empty if pi == key
                if (pi != key)
                {
                    // scan A(i,:)
                    for (GrB_Index pS = AP (i) ; pS < AP (i+1) ; pS++)
                    {
                        // get A(i,j)
                        int64_t j = AJ (pS) ;
                        if (PARENT (j) != key)
                        {
                            // add the entry T(i,j) to T, but skip it if
                            // PARENT (j) is equal to key
                            SET_TJ (p, j)       // Tj [p] = j ;
                            p++ ;
                        }
                    }
                    // Add the entry T(i,key) if there is room for it in T(i,:);
                    // if and only if node i is adjacent to a node j in the
                    // largest component.  The only way there can be space if
                    // at least one T(i,j) appears with PARENT (j) equal to the
                    // key (that is, node j is in the largest connected
                    // component, key == PARENT (j).  One of these j's can then
                    // be replaced with the key.  If node i is not adjacent to
                    // any node in the largest component, then there is no
                    // space in T(i,:) and no new edge to the largest component
                    // is added.
                    if (p - pstart < AP (i+1) - AP (i))
                    {
                        SET_TJ (p, key) ;       // Tj [p] = key ;
                        p++ ;
                    }
                }
            }
            // count the number of entries inserted into T by this thread
            count [tid] = p - TP (ktid) ;
        }

        // Compact empty space out of Tj not filled in from the above phase.
        nvals = 0 ;
        for (tid = 0 ; tid < nthreads ; tid++)
        {
            int64_t ktid = range [tid]  ;
            memmove (Tj32 ? ((void *) (Tj32 + nvals)) : ((void *) (Tj64 + nvals)),
                     Tj32 ? ((void *) (Tj32 + TP (ktid))) : ((void *) (Tj64 + TP (ktid))),
                     tjsize * count [tid]) ;

#if 0
            if (Tj32)
            {
                memmove (Tj32 + nvals, Tj32 + TP (ktid),
                    sizeof (uint32_t) * count [tid]) ;
            }
            else
            {
                memmove (Tj64 + nvals, Tj64 + TP (ktid),
                    sizeof (uint64_t) * count [tid]) ;
            }
#endif

            nvals += count [tid] ;
            count [tid] = nvals - count [tid] ;
        }

        // Compact empty space out of Tp
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (tid = 0 ; tid < nthreads ; tid++)
        {
            int64_t ktid = range [tid] ;
            GrB_Index p = TP (ktid) ;
            for (int64_t i = range [tid] ; i < range [tid+1] ; i++)
            {
                int64_t tp = TP (i) ;
                tp -= p - count [tid] ;
                SET_TP (i, tp) ;            // Tp [i] = tp ;
            }
        }

        // finalize T
        SET_TP (n, nvals) ;     // Tp [n] = nvals ;
        Tj_len = nvals ;

        // load T_Container->p,i from the C arrays, Tp and Tj, for final phase
        GRB_TRY (GxB_Vector_load (T_Container->p, &Tp, Tp_type, Tp_len,
            Tp_size, GrB_DEFAULT, NULL)) ;
        GRB_TRY (GxB_Vector_load (T_Container->i, &Tj, Tj_type, Tj_len,
            Tj_size, GrB_DEFAULT, NULL)) ;

        T_Container->nrows_nonempty = -1 ;
        T_Container->ncols_nonempty = -1 ;
        T_Container->jumbled = true ;
        T_Container->nvals = nvals ;

        // load T in sparse CSR format from the T_Container
        GRB_TRY (GxB_load_Matrix_from_Container (T, T_Container, NULL)) ;

        // load A_Container->p,i from the C arrays, Ap and Aj
        // This is the original G->A, and it is unchanged.
        GRB_TRY (GxB_Vector_load (A_Container->p, &Ap, Ap_type, Ap_len,
            Ap_size, Ap_handling, NULL)) ;
        GRB_TRY (GxB_Vector_load (A_Container->i, &Aj, Aj_type, Aj_len,
            Aj_size, Aj_handling, NULL)) ;

        // load A in sparse CSR format from the A_Container
        GRB_TRY (GxB_load_Matrix_from_Container (A, A_Container, NULL)) ;

        //----------------------------------------------------------------------
        // load the Px array back into the parent vector
        //----------------------------------------------------------------------

        GRB_TRY (GxB_Vector_load (parent, &Px, type, n, Px_size,
            GrB_DEFAULT, NULL)) ;

// ].  The unload/load of A into Ap, Aj, Ax will not be needed, and G->A
// will become truly a read-only matrix.

        // final phase uses the pruned matrix T
        A = T ;

        #ifdef TIMINGS
        toc2 = LAGraph_WallClockTime ( ) ;
        timings [7] = toc2 - tic2 ;  // prune
        tic2 = toc2 ;
        #endif
    }

    #ifdef TIMINGS
    toc = LAGraph_WallClockTime ( ) ;
    timings [1] = toc - tic ;  // total sampling time
    tic = toc ;
    #endif

    //--------------------------------------------------------------------------
    // check for quick return
    //--------------------------------------------------------------------------

    // The sample phase may have already found that G->A has a single component,
    // in which case the matrix A is now empty.

    if (nvals == 0)
    {
        (*component) = parent ;
        LG_FREE_WORK ;
        #ifdef TIMINGS
        print_timings (timings) ;
        LG_SET_BURBLE (false) ;
        #endif
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // final phase
    //--------------------------------------------------------------------------

    LG_TRY (fastsv (A, parent2, mngp, &gp, &gp_new, t, eq, min, min_2nd,
        parent, ramp_v, msg)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*component) = parent ;
    parent = NULL ;
    LG_FREE_WORK ;
    #ifdef TIMINGS
    toc = LAGraph_WallClockTime ( ) ;
    timings [2] = toc - tic ;  // final phase
    print_timings (timings) ;
    LG_SET_BURBLE (false) ;
    #endif
    return (GrB_SUCCESS) ;
#else
    LG_ASSERT (false, GrB_NOT_IMPLEMENTED) ;
#endif
}
