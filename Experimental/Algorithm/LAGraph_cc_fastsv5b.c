//------------------------------------------------------------------------------
// LAGraph_cc_fastsv5b: connected components
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2020 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

/**
 * Code is based on the algorithm described in the following paper
 * Zhang, Azad, Hu. FastSV: FastSV: A Distributed-Memory Connected Component
 * Algorithm with Fast Convergence (SIAM PP20)
 *
 * Modified by Tim Davis, Texas A&M University
 **/

// The input matrix A must be symmetric.  Self-edges (diagonal entries) are
// OK, and are ignored.  The values and type of A are ignored; just its
// pattern is accessed.

// The matrix A must have dimension 2^32 or less.  If it is larger, use the
// 64-bit version of this method instead.  TODO combine the two versions into a
// single user-callable code.

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// Reduce_assign32:  w (index) += src, using MIN as the "+=" accum operator
//------------------------------------------------------------------------------

// mask = NULL, accumulator = GrB_MIN_UINT32, descriptor = NULL.
// Duplicates are summed with the accumulator, which differs from how
// GrB_assign works.  GrB_assign states that the presence of duplicates results
// in undefined behavior.  SuiteSparse:GraphBLAS follows the MATLAB rule, which
// discards all but the first of the duplicates.  TODO: add this to GraphBLAS
// as a variant of GrB_assign, either as GxB_assign_accum (or another name),
// or as a GxB_* descriptor setting.

#define LAGRAPH_FREE_ALL

// hash table
const int P = 1024;

#define HASH(x) (((x << 4) + x) & 1023)
#define NEXT(x) ((x + 23) & 1023)

static inline void ht_init (int *ht_key, int *ht_val)
{
    memset(ht_key, -1, sizeof(int) * P);
    memset(ht_val,  0, sizeof(int) * P);
}

static inline void ht_sample (uint32_t *V32, int n, int samples, int *ht_key, int *ht_val)
{
    for (int i = 0; i < samples; i++) {
        int x = V32 [rand() % n];
        int h = HASH (x);
        while (ht_key [h] != -1 && ht_key [h] != x)
            h = NEXT (h);
        ht_key [h] = x;
        ht_val [h] += 1;
    }
}

static inline int ht_most_frequent (int *ht_key, int *ht_val)
{
    int key = -1, val = 0;
    for (int i = 0; i < P; i++)
        if (ht_val [i] > val)
        {
            key = ht_key [i];
            val = ht_val [i];
        }
    return key;
}

static inline GrB_Info Reduce_assign32
(
    GrB_Vector *w_handle,   // vector of size n, all entries present
    GrB_Vector *s_handle,   // vector of size n, all entries present
    uint32_t *index,        // array of size n
    GrB_Index n,
    int nthreads,
    int *ht_key,
    int *ht_val
)
{
    GrB_Type w_type, s_type ;
    GrB_Index w_n, s_n, w_nvals, s_nvals, *w_i, *s_i, w_size, s_size ;
    uint32_t *w_x, *s_x ;

    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,1)
    LAGr_Vector_export_Full (w_handle, &w_type, &w_n, (void **) &w_x, 
        &w_size, NULL) ;
    LAGr_Vector_export_Full (s_handle, &s_type, &s_n, (void **) &s_x,
        &s_size, NULL) ;
    #elif GxB_IMPLEMENTATION == GxB_VERSION (4,0,0)
    LAGr_Vector_export_Full (w_handle, &w_type, &w_n, (void **) &w_x, NULL) ;
    LAGr_Vector_export_Full (s_handle, &s_type, &s_n, (void **) &s_x, NULL) ;
    #else
    LAGr_Vector_export (w_handle, &w_type, &w_n, &w_nvals, &w_i,
        (void **) &w_x, NULL) ;
    LAGr_Vector_export (s_handle, &s_type, &s_n, &s_nvals, &s_i,
        (void **) &s_x, NULL) ;
    #endif

    if (nthreads >= 4)
    {
        uint32_t *mem = LAGraph_malloc (nthreads * P, sizeof (uint32_t));

        ht_init (ht_key, ht_val) ;
        ht_sample (index, n, 864, ht_key, ht_val) ;

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int t = 0; t < nthreads; t++)
        {
            uint32_t *buf = mem + t * P;
            for (int h = 0; h < P; h++)
            {
                if (ht_key [h] != -1)
                {
                    buf [h] = w_x [ht_key [h]];
                }
            }
            int st = (n * t + nthreads - 1) / nthreads;
            int ed = (n * t + n + nthreads - 1) / nthreads;

            for (int k = st ; k < ed ; k++)
            {
                uint32_t i = index [k] ;
                int h = HASH(i);
                while (ht_key [h] != -1 && ht_key [h] != i)
                {
                    h = NEXT (h);
                }

                if (ht_key [h] == -1)
                {
                    w_x [i] = LAGRAPH_MIN (w_x [i], s_x [k]);
                }
                else
                {
                    buf [h] = LAGRAPH_MIN (buf [h], s_x [k]);
                }
            }
        }

        for (int h = 0; h < P; h++)
        {
            int i = ht_key [h];
            if (i != -1)
            {
                for (int j = 0; j < nthreads; j++)
                {
                    w_x [i] = LAGRAPH_MIN (w_x [i], mem [j * P + h]);
                }
            }
        }

        LAGRAPH_FREE (mem);
    }
    else
    {
        // sequential version, to avoid atomics
        for (GrB_Index k = 0 ; k < n ; k++)
        {
            uint32_t i = index [k] ;
            w_x [i] = LAGRAPH_MIN (w_x [i], s_x [k]) ;
        }
    }

    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,1)
    LAGr_Vector_import_Full (w_handle, w_type, w_n, (void **) &w_x, 
        w_size, NULL) ;
    LAGr_Vector_import_Full (s_handle, s_type, s_n, (void **) &s_x,
        s_size, NULL) ;
    #elif GxB_IMPLEMENTATION == GxB_VERSION (4,0,0)
    LAGr_Vector_import_Full (w_handle, w_type, w_n, (void **) &w_x, NULL) ;
    LAGr_Vector_import_Full (s_handle, s_type, s_n, (void **) &s_x, NULL) ;
    #else
    LAGr_Vector_import (w_handle, w_type, w_n, w_nvals, &w_i,
        (void **) &w_x, NULL) ;
    LAGr_Vector_import (s_handle, s_type, s_n, s_nvals, &s_i,
        (void **) &s_x, NULL) ;
    #endif
    return (GrB_SUCCESS) ;
}

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE (I) ;              \
    LAGRAPH_FREE (V32) ;            \
    LAGRAPH_FREE (ht_key) ;         \
    LAGRAPH_FREE (ht_val) ;         \
    LAGr_free (&f) ;                \
    LAGr_free (&gp) ;               \
    LAGr_free (&mngp) ;             \
    LAGr_free (&gp_new) ;           \
    LAGr_free (&mod) ;              \
}

//------------------------------------------------------------------------------
// LAGraph_cc_fastsv5
//------------------------------------------------------------------------------

GrB_Info LAGraph_cc_fastsv5b
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix *A,          // input matrix
                            //   content remains the same, but pointer changes
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info ;
    uint32_t *V32 = NULL ;
    int *ht_key = NULL, *ht_val = NULL;
    GrB_Index n, nnz, *I = NULL ;
    GrB_Vector f = NULL, gp_new = NULL, mngp = NULL, mod = NULL, gp = NULL ;
    GrB_Matrix S = NULL, T = NULL ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LAGr_Matrix_nrows (&n, *A) ;
    LAGr_Matrix_nvals (&nnz, *A) ;

    if (n > UINT32_MAX)
    {
        LAGRAPH_ERROR ("problem too large; use 64-bit version instead",
            GrB_INVALID_VALUE) ;
    }

#define FASTSV_SAMPLES 4

    GxB_Format_Value format;
    LAGRAPH_OK (GxB_get (*A , GxB_FORMAT, &format)) ;
    bool sampling = (format == GxB_BY_ROW) && (n * FASTSV_SAMPLES * 2 < nnz);
    
    if (sanitize)
    {
        // S = A | A'
        LAGr_Matrix_new (&S, GrB_BOOL, n, n) ;
        LAGr_eWiseAdd (S, NULL, NULL, GrB_LOR, *A, *A, LAGraph_desc_otoo) ;
    }
    else
    {
        // Use the input as-is, and assume it is symmetric
        // LAGr_Matrix_dup (&S, A) ;
        S = *A;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // determine # of threads to use for Reduce_assign
    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (nthreads, n / 16) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    // # of threads to use for typecast
    int nthreads2 = n / (64*1024) ;
    nthreads2 = LAGRAPH_MIN (nthreads2, nthreads) ;
    nthreads2 = LAGRAPH_MAX (nthreads2, 1) ;
    // printf ("nthreads %d nthreads2 %d\n", nthreads, nthreads2) ;

    // vectors
    LAGr_Vector_new (&f,      GrB_UINT32, n) ;
    LAGr_Vector_new (&gp_new, GrB_UINT32, n) ;
    LAGr_Vector_new (&mod,    GrB_BOOL, n) ;
    // temporary arrays
    I   = LAGraph_malloc (n, sizeof (GrB_Index)) ;
    V32 = LAGraph_malloc (n, sizeof (uint32_t)) ;
    // prepare vectors
    #pragma omp parallel for num_threads(nthreads2) schedule(static)
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        I [i] = i ;
        V32 [i] = (uint32_t) i ;
    }
    LAGr_Vector_build (f, I, V32, n, GrB_PLUS_UINT32) ;
    LAGr_Vector_dup (&gp,   f) ;
    LAGr_Vector_dup (&mngp, f) ;

    ht_key = LAGraph_malloc (P, sizeof (int));
    ht_val = LAGraph_malloc (P, sizeof (int));

    //--------------------------------------------------------------------------
    // main computation
    //--------------------------------------------------------------------------

    if (sampling)
    {
        GrB_Type type;
        GrB_Index nrows, ncols, nvals;
        int64_t nonempty;
        GrB_Index *Sp, *Sj;
        void *Sx;
        bool S_jumbled = false ;
        GrB_Index Sp_size, Sj_size, Sx_size ;

        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,1)
        GrB_Matrix_nvals (&nvals, S) ;
        GxB_Matrix_export_CSR (&S, &type, &nrows, &ncols,
                &Sp, &Sj, &Sx, &Sp_size, &Sj_size, &Sx_size, &S_jumbled, NULL) ;
        #elif GxB_IMPLEMENTATION == GxB_VERSION (4,0,0)
        GxB_Matrix_export_CSR (&S, &type, &nrows, &ncols, &nvals,
                &S_jumbled, &nonempty, &Sp, &Sj, &Sx, NULL);
        #else
        GxB_Matrix_export_CSR (&S, &type, &nrows, &ncols, &nvals,
                &nonempty, &Sp, &Sj, &Sx, NULL);
        #endif

        GrB_Index Tp_size = nrows+1 ;
        GrB_Index Tj_size = nvals ;
        GrB_Index Tx_size = nvals ;

        GrB_Index *Tp = LAGraph_malloc (Tp_size, sizeof (GrB_Index)) ;
        GrB_Index *Tj = LAGraph_malloc (Tj_size, sizeof (GrB_Index)) ;
        void *Tx = LAGraph_malloc (Tx_size, 1) ;

        int *range = LAGraph_malloc (nthreads + 1, sizeof (int)) ;
        GrB_Index *count = LAGraph_malloc (nthreads + 1, sizeof (GrB_Index)) ;
        memset (count,  0, sizeof (GrB_Index) * (nthreads + 1)) ;

        for (int i = 0; i <= nthreads; i++)
        {
            range [i] = (n * i + nthreads - 1) / nthreads;
        }

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int t = 0; t < nthreads; t++)
        {
            for (int i = range[t]; i < range[t + 1]; i++)
            {
                int deg = Sp [i + 1] - Sp [i];
                count [t + 1] += LAGRAPH_MIN (FASTSV_SAMPLES, deg) ;
            }
        }

        for (int i = 0; i < nthreads; i++)
        {
            count [i + 1] += count [i];
        }

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int t = 0; t < nthreads; t++)
        {
            GrB_Index p = count [t];
            Tp [range [t]] = p;
            for (int i = range[t]; i < range[t + 1]; i++)
            {
                for (int j = 0; j < FASTSV_SAMPLES && Sp [i] + j < Sp [i + 1];
                    j++)
                {
                    Tj [p++] = Sj [Sp [i] + j];
                }
                Tp [i + 1] = p;
            }
        }

        GrB_Index t_nvals = Tp[nrows];

        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,1)
        GxB_Matrix_import_CSR (&T, type, nrows, ncols,
                &Tp, &Tj, &Tx, Tp_size, Tj_size, Tx_size, S_jumbled, NULL) ;
        #elif GxB_IMPLEMENTATION == GxB_VERSION (4,0,0)
        GxB_Matrix_import_CSR (&T, type, nrows, ncols, t_nvals,
                S_jumbled, -1, &Tp, &Tj, &Tx, NULL);
        #else
        GxB_Matrix_import_CSR (&T, type, nrows, ncols, t_nvals,
                -1, &Tp, &Tj, &Tx, NULL);
        #endif

        bool change = true, is_first = true;
        while (change)
        {
            // hooking & shortcutting
            LAGr_mxv (mngp, NULL, GrB_MIN_UINT32, GxB_MIN_SECOND_UINT32, T, gp,
                NULL) ;
            if (!is_first)
            {
                LAGRAPH_OK (Reduce_assign32 (&f, &mngp, V32, n, nthreads, ht_key, ht_val)) ;
            }
            // old:
            // LAGr_eWiseMult (f, NULL, NULL, GrB_MIN_UINT32, f, mngp, NULL) ;
            // LAGr_eWiseMult (f, NULL, NULL, GrB_MIN_UINT32, f, gp, NULL) ;
            // new:
            LAGr_eWiseAdd (f, NULL, GrB_MIN_UINT32, GrB_MIN_UINT32, mngp, gp,
                NULL) ;
            // calculate grandparent
            LAGr_Vector_extractTuples (NULL, V32, &n, f) ;
            #pragma omp parallel for num_threads(nthreads2) schedule(static)
            for (uint32_t i = 0 ; i < n ; i++)
            {
                I [i] = (GrB_Index) V32 [i] ;
            }
            LAGr_extract (gp_new, NULL, NULL, f, I, n, NULL) ;
            // check termination
            LAGr_eWiseMult (mod, NULL, NULL, GrB_NE_UINT32, gp_new, gp, NULL) ;
            LAGr_reduce (&change, NULL, GxB_LOR_BOOL_MONOID, mod, NULL) ;
            // swap gp and gp_new
            GrB_Vector t = gp ; gp = gp_new ; gp_new = t ;
            is_first = false;
        }

        ht_init(ht_key, ht_val) ;
        ht_sample (V32, n, 864, ht_key, ht_val) ;
        int key = ht_most_frequent(ht_key, ht_val) ;

        int64_t t_nonempty = -1 ;
        bool T_jumbled = false ;

        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,1)
        GxB_Matrix_export_CSR (&T, &type, &nrows, &ncols,
                &Tp, &Tj, &Tx, &Tp_size, &Tj_size, &Tx_size,
                &T_jumbled, NULL) ;
        #elif GxB_IMPLEMENTATION == GxB_VERSION (4,0,0)
        GxB_Matrix_export_CSR (&T, &type, &nrows, &ncols, &t_nvals,
                &T_jumbled, &t_nonempty, &Tp, &Tj, &Tx, NULL);
        #else
        GxB_Matrix_export_CSR (&T, &type, &nrows, &ncols, &t_nvals,
                &t_nonempty, &Tp, &Tj, &Tx, NULL);
        #endif

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int t = 0; t < nthreads; t++)
        {
            GrB_Index ptr = Sp[range[t]];
            for (int v = range[t]; v < range[t + 1]; v++)
            {
                int pv = V32 [v];
                Tp [v] = ptr;
                if (pv != key)
                {
                    for (GrB_Index i = Sp [v]; i < Sp [v + 1]; i++)
                    {
                        int u = Sj [i];
                        if (V32 [u] != key)
                        {
                            Tj [ptr++] = u;
                        }
                    }
                    if (ptr - Tp[v] < Sp [v + 1] - Sp [v])
                    {
                        Tj [ptr++] = key;
                    }
                }
            }
            count[t] = ptr - Tp [range [t]];
        }

        GrB_Index offset = 0;
        for (int i = 0; i < nthreads; i++)
        {
            memcpy(Tj + offset, Tj + Tp [range [i]], sizeof(GrB_Index) * count[i]);
            offset += count[i];
            count[i] = offset - count[i];
        }

        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int t = 0; t < nthreads; t++)
        {
            GrB_Index ptr = Tp [range [t]];
            for (int v = range[t]; v < range[t + 1]; v++)
            {
                Tp [v] -= ptr - count[t];
            }
        }
        Tp [n] = offset;
        LAGRAPH_FREE (count);
        LAGRAPH_FREE (range);

        #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,1)
        GxB_Matrix_import_CSR (&S, type, nrows, ncols,
                &Sp, &Sj, &Sx, Sp_size, Sj_size, Sx_size, S_jumbled, NULL);
        GxB_Matrix_import_CSR (&T, type, nrows, ncols,
                &Tp, &Tj, &Tx, Tp_size, Tj_size, Tx_size, T_jumbled, NULL) ;
        #elif GxB_IMPLEMENTATION == GxB_VERSION (4,0,0)
        GxB_Matrix_import_CSR (&S, type, nrows, ncols, nvals,
                S_jumbled, nonempty, &Sp, &Sj, &Sx, NULL);
        GxB_Matrix_import_CSR (&T, type, nrows, ncols, offset,
                T_jumbled, -1, &Tp, &Tj, &Tx, NULL);
        #else
        GxB_Matrix_import_CSR (&S, type, nrows, ncols, nvals,
                nonempty, &Sp, &Sj, &Sx, NULL);
        GxB_Matrix_import_CSR (&T, type, nrows, ncols, offset,
                -1, &Tp, &Tj, &Tx, NULL);
        #endif

    }
    else
    {
        T = S;
    }

    LAGr_Matrix_nvals (&nnz, T);
    bool change = true;
    while (change && nnz > 0)
    {
        // hooking & shortcutting
        LAGr_mxv (mngp, NULL, GrB_MIN_UINT32, GxB_MIN_SECOND_UINT32, T, gp,
            NULL) ;
        LAGRAPH_OK (Reduce_assign32 (&f, &mngp, V32, n, nthreads, ht_key, ht_val)) ;
        // old:
        // LAGr_eWiseMult (f, NULL, NULL, GrB_MIN_UINT32, f, mngp, NULL) ;
        // LAGr_eWiseMult (f, NULL, NULL, GrB_MIN_UINT32, f, gp, NULL) ;
        // new:
        LAGr_eWiseAdd (f, NULL, GrB_MIN_UINT32, GrB_MIN_UINT32, mngp, gp, NULL);
        // calculate grandparent

        LAGr_Vector_extractTuples (NULL, V32, &n, f) ;
        #pragma omp parallel for num_threads(nthreads2) schedule(static)
        for (uint32_t i = 0 ; i < n ; i++)
        {
            I [i] = (GrB_Index) V32 [i] ;
        }
        LAGr_extract (gp_new, NULL, NULL, f, I, n, NULL) ;
        // check termination
        LAGr_eWiseMult (mod, NULL, NULL, GrB_NE_UINT32, gp_new, gp, NULL) ;
        LAGr_reduce (&change, NULL, GxB_LOR_BOOL_MONOID, mod, NULL) ;
        // swap gp and gp_new
        GrB_Vector t = gp ; gp = gp_new ; gp_new = t ; 
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    *result = f ;
    f = NULL ;
    if (!sanitize)
    {
        *A = S;
    }
    else
    {
        LAGr_free (&S) ;
    }
    if (sampling)
    {
        LAGr_free (&T) ;
    }
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

