//------------------------------------------------------------------------------
// LAGraph_cc_fastsv4: connected components
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

#include "LAGraph.h"

//------------------------------------------------------------------------------
// atomic_min_uint32: compute (*p) = min (*p, value), via atomic update
//------------------------------------------------------------------------------

static inline void atomic_min_uint32
(
    uint32_t *p,        // input/output
    uint32_t value      // input
)
{
    uint32_t old, new ;
    do
    {
        // get the old value at (*p)
        // #pragma omp atomic read
        old = (*p) ;
        // compute the new minimum
        new = LAGRAPH_MIN (old, value) ;
    }
    while (!__sync_bool_compare_and_swap (p, old, new)) ;
}

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

static GrB_Info Reduce_assign32
(
    GrB_Vector *w_handle,   // vector of size n, all entries present
    GrB_Vector *s_handle,   // vector of size n, all entries present
    uint32_t *index,        // array of size n
    GrB_Index n,
    int nthreads
)
{

    GrB_Type w_type, s_type ;
    GrB_Index w_n, s_n, w_nvals, s_nvals, *w_i, *s_i ;
    uint32_t *w_x, *s_x ;

    LAGr_Vector_export (w_handle, &w_type, &w_n, &w_nvals, &w_i,
        (void **) &w_x, NULL) ;
    LAGr_Vector_export (s_handle, &s_type, &s_n, &s_nvals, &s_i,
        (void **) &s_x, NULL) ;

    if (nthreads >= 4)
    {
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (GrB_Index k = 0 ; k < n ; k++)
        {
            uint32_t i = index [k] ;
            atomic_min_uint32 (&(w_x [i]), s_x [k]) ;
        }
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

    LAGr_Vector_import (w_handle, w_type, w_n, w_nvals, &w_i,
        (void **) &w_x, NULL) ;
    LAGr_Vector_import (s_handle, s_type, s_n, s_nvals, &s_i,
        (void **) &s_x, NULL) ;

    return (GrB_SUCCESS) ;
}

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE (I) ;              \
    LAGRAPH_FREE (V32) ;            \
    LAGr_free (&f) ;                \
    LAGr_free (&gp) ;               \
    LAGr_free (&mngp) ;             \
    LAGr_free (&gp_new) ;           \
    LAGr_free (&mod) ;              \
    if (sanitize) LAGr_free (&S) ;  \
}

//------------------------------------------------------------------------------
// LAGraph_cc_fastsv4
//------------------------------------------------------------------------------

GrB_Info LAGraph_cc_fastsv4
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info ;
    uint32_t *V32 = NULL ;
    GrB_Index n, *I = NULL ;
    GrB_Vector f = NULL, gp_new = NULL, mngp = NULL, mod = NULL, gp = NULL ;
    GrB_Matrix S = NULL ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LAGr_Matrix_nrows (&n, A) ;

    if (n > UINT32_MAX)
    {
        LAGRAPH_ERROR ("problem too large; use 64-bit version instead",
            GrB_INVALID_VALUE) ;
    }

    if (sanitize)
    {
        // S = A | A'
        LAGr_Matrix_new (&S, GrB_BOOL, n, n) ;
        LAGr_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, LAGraph_desc_otoo) ;
    }
    else
    {
        // Use the input as-is, and assume it is symmetric
        S = A ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // determine # of threads to use for Reduce_assign
    int nthreads_max = LAGraph_get_nthreads ( ) ;
    int nthreads = n / (1024*1024) ;
    nthreads = LAGRAPH_MIN (nthreads, nthreads_max) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    // # of threads to use for typecast
    int nthreads2 = n / (64*1024) ;
    nthreads2 = LAGRAPH_MIN (nthreads2, nthreads_max) ;
    nthreads2 = LAGRAPH_MAX (nthreads2, 1) ;

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

    //--------------------------------------------------------------------------
    // main computation
    //--------------------------------------------------------------------------

    bool diff = true ;
    while (diff)
    {
        // hooking & shortcutting
        LAGr_mxv (mngp, NULL, GrB_MIN_UINT32, GxB_MIN_SECOND_UINT32, S, gp,
            NULL) ;
        LAGRAPH_OK (Reduce_assign32 (&f, &mngp, V32, n, nthreads)) ;
        LAGr_eWiseMult (f, NULL, NULL, GrB_MIN_UINT32, f, mngp, NULL) ;
        LAGr_eWiseMult (f, NULL, NULL, GrB_MIN_UINT32, f, gp, NULL) ;
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
        LAGr_reduce (&diff, NULL, GxB_LOR_BOOL_MONOID, mod, NULL) ;
        // swap gp and gp_new
        GrB_Vector t = gp ; gp = gp_new ; gp_new = t ; 
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    *result = f ;
    f = NULL ;
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

