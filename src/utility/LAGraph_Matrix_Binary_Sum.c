//------------------------------------------------------------------------------
// LAGraph_Matrix_Binary_Sum: sum an array of matrices by binary reduction
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Michel Pelletier.

//------------------------------------------------------------------------------

// LAGraph_Matrix_Binary_Sum combines an array of matrices into a single matrix
// C, computing the same result as LAGraph_Matrix_Sum but with a different
// technique: a pairwise binary reduction tree built from GrB_eWiseAdd, as
// described in the GraphChallenge "Anonymized Network Sensing" paper (Fig. 3,
// "Binary Summation of Traffic Matrices").  At each level of the tree the
// matrices are summed in disjoint adjacent pairs; an unpaired (odd) trailing
// matrix is carried up to the next level unchanged.  The levels repeat until a
// single matrix remains.  Each pair-sum uses the binary operator dup, which has
// the same set-union semantics as the dup operator of GrB_Matrix_build: dup is
// applied wherever both inputs have an entry, and lone entries pass through.
// With dup = GrB_PLUS_FP64 (for example) this computes the element-wise sum of
// all the input matrices.

// The independent pair-sums within a level are issued concurrently across
// LG_nthreads_outer threads with OpenMP; each such outer thread calls into
// GraphBLAS, which parallelizes each GrB_eWiseAdd internally with
// LG_nthreads_inner threads (the documented two-level model).

// Memory: a level's intermediate results are only needed as inputs to the next
// level, so each pair's two inputs are freed as soon as their sum is built.
// Only the active merge frontier (the current level plus the partially-built
// next level) is ever resident, rather than every intermediate result at once.
// This is what makes the tree a memory win: adding two matrices each with N
// entries yields a matrix with fewer than 2N entries, and the smaller
// intermediates stay in faster memory.  The original input matrices (level 0)
// are never freed.

// Ownership is tracked per working-array slot: ownedW[k] is true iff W[k] is a
// matrix this function created (and so must free), and false for an original
// input.  This makes LG_FREE_WORK a simple sweep that frees exactly the
// matrices we still own, and is correct at any point -- including partial
// failure inside a level, since GrB_free of an already-freed (NULL) handle is a
// no-op.  Freeing a consumed input sets its handle to NULL, so the same slot is
// never freed twice.

// All input matrices must have identical dimensions and identical built-in
// type; C is created with that same type and dimensions.  Unlike
// LAGraph_Matrix_Sum, the dup operator must be non-NULL, since GrB_eWiseAdd
// requires a binary operator.

#define LG_FREE_WORK                                        \
{                                                           \
    if (W != NULL && ownedW != NULL)                        \
    {                                                       \
        for (GrB_Index k = 0 ; k < nmatrices ; k++)         \
        {                                                   \
            if (ownedW [k]) GrB_free (& W [k]) ;            \
        }                                                   \
    }                                                       \
    if (Wnext != NULL && ownedN != NULL)                    \
    {                                                       \
        for (GrB_Index k = 0 ; k < nmatrices ; k++)         \
        {                                                   \
            if (ownedN [k]) GrB_free (& Wnext [k]) ;        \
        }                                                   \
    }                                                       \
    LAGraph_Free ((void **) &W, NULL) ;                     \
    LAGraph_Free ((void **) &ownedW, NULL) ;                \
    LAGraph_Free ((void **) &Wnext, NULL) ;                 \
    LAGraph_Free ((void **) &ownedN, NULL) ;                \
}

#define LG_FREE_ALL                                         \
{                                                           \
    LG_FREE_WORK ;                                          \
    GrB_free (C) ;                                          \
}

#include "LG_internal.h"

int LAGraph_Matrix_Binary_Sum
(
    // output:
    GrB_Matrix *C,          // result = combination of all input matrices
    // input:
    GrB_Matrix *Matrices,   // array of nmatrices input matrices
    GrB_Index nmatrices,    // number of matrices in the array (must be >= 1)
    GrB_BinaryOp dup,       // operator to combine (i,j) entries (must be != NULL)
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix *W = NULL, *Wnext = NULL ;
    bool *ownedW = NULL, *ownedN = NULL ;
    LG_ASSERT_MSG (C != NULL, GrB_NULL_POINTER, "&C != NULL") ;
    LG_ASSERT (Matrices != NULL, GrB_NULL_POINTER) ;
    (*C) = NULL ;
    LG_ASSERT_MSG (nmatrices >= 1, GrB_INVALID_VALUE,
        "nmatrices must be >= 1") ;
    LG_ASSERT_MSG (dup != NULL, GrB_NULL_POINTER,
        "dup operator must be non-NULL") ;
    LG_ASSERT (Matrices [0] != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the reference dimensions and type from the first matrix
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    int32_t typecode ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, Matrices [0])) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, Matrices [0])) ;
    GRB_TRY (GrB_get (Matrices [0], &typecode, GrB_EL_TYPE_CODE)) ;

    // map the type code to a built-in GrB_Type for the intermediate matrices;
    // GrB_eWiseAdd itself is polymorphic, so only GrB_Matrix_new needs the type
    GrB_Type gtype = NULL ;
    switch (typecode)
    {
        case GrB_BOOL_CODE   : gtype = GrB_BOOL   ; break ;
        case GrB_INT8_CODE   : gtype = GrB_INT8   ; break ;
        case GrB_INT16_CODE  : gtype = GrB_INT16  ; break ;
        case GrB_INT32_CODE  : gtype = GrB_INT32  ; break ;
        case GrB_INT64_CODE  : gtype = GrB_INT64  ; break ;
        case GrB_UINT8_CODE  : gtype = GrB_UINT8  ; break ;
        case GrB_UINT16_CODE : gtype = GrB_UINT16 ; break ;
        case GrB_UINT32_CODE : gtype = GrB_UINT32 ; break ;
        case GrB_UINT64_CODE : gtype = GrB_UINT64 ; break ;
        case GrB_FP32_CODE   : gtype = GrB_FP32   ; break ;
        case GrB_FP64_CODE   : gtype = GrB_FP64   ; break ;
        default :
            LG_ASSERT_MSG (false, GrB_NOT_IMPLEMENTED,
                "user-defined types are not supported") ;
    }

    //--------------------------------------------------------------------------
    // validate every matrix has the same dimensions and type
    //--------------------------------------------------------------------------

    for (GrB_Index k = 0 ; k < nmatrices ; k++)
    {
        GrB_Matrix Ak = Matrices [k] ;
        LG_ASSERT (Ak != NULL, GrB_NULL_POINTER) ;
        GrB_Index r, c ;
        int32_t code ;
        GRB_TRY (GrB_Matrix_nrows (&r, Ak)) ;
        GRB_TRY (GrB_Matrix_ncols (&c, Ak)) ;
        LG_ASSERT_MSG (r == nrows && c == ncols, GrB_DIMENSION_MISMATCH,
            "all input matrices must have the same dimensions") ;
        GRB_TRY (GrB_get (Ak, &code, GrB_EL_TYPE_CODE)) ;
        LG_ASSERT_MSG (code == typecode, GrB_DOMAIN_MISMATCH,
            "all input matrices must have the same type") ;
    }

    //--------------------------------------------------------------------------
    // handle the single-matrix case: C is an independent copy of Matrices [0]
    //--------------------------------------------------------------------------

    if (nmatrices == 1)
    {
        GRB_TRY (GrB_Matrix_dup (C, Matrices [0])) ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // allocate the double-buffered working arrays and their ownership flags
    //--------------------------------------------------------------------------

    // W / ownedW hold the current level; Wnext / ownedN hold the level being
    // built.  ownedW is initialized (to all-false, level 0 being the original
    // inputs) before Wnext is allocated, so that if a later allocation fails the
    // LG_FREE_WORK sweep never reads an uninitialized flag array.

    LG_TRY (LAGraph_Malloc ((void **) &W, nmatrices, sizeof (GrB_Matrix),
        msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &ownedW, nmatrices, sizeof (bool),
        msg)) ;
    for (GrB_Index k = 0 ; k < nmatrices ; k++)
    {
        W [k] = Matrices [k] ;      // level 0 references the original inputs
        ownedW [k] = false ;        // which we never free
    }
    LG_TRY (LAGraph_Malloc ((void **) &Wnext, nmatrices, sizeof (GrB_Matrix),
        msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &ownedN, nmatrices, sizeof (bool),
        msg)) ;
    for (GrB_Index k = 0 ; k < nmatrices ; k++)
    {
        Wnext [k] = NULL ;
        ownedN [k] = false ;
    }

    //--------------------------------------------------------------------------
    // pairwise binary reduction
    //--------------------------------------------------------------------------

    GrB_Index count = nmatrices ;   // number of matrices at the current level

    while (count > 1)
    {
        GrB_Index npairs = count / 2 ;          // number of pair-sums
        int64_t p ;

        // reset the next-level buffers for this level
        for (GrB_Index k = 0 ; k < nmatrices ; k++)
        {
            Wnext [k] = NULL ;
            ownedN [k] = false ;
        }

        // outer threads for the independent pair-sums at this level; each calls
        // GraphBLAS, which nests inner threads underneath (two-level model)
        int nthreads = LG_nthreads_outer ;
        nthreads = LAGRAPH_MIN (nthreads, (int) npairs) ;
        nthreads = LAGRAPH_MAX (nthreads, 1) ;

        // GRB_TRY cannot be used inside the OpenMP region (it returns from the
        // function), so the first error is captured into sum_status under a
        // critical section and checked after the region.
        int sum_status = GrB_SUCCESS ;

        #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
        for (p = 0 ; p < (int64_t) npairs ; p++)
        {
            GrB_Matrix R = NULL ;
            GrB_Info info = GrB_Matrix_new (&R, gtype, nrows, ncols) ;
            if (info >= GrB_SUCCESS)
            {
                info = GrB_eWiseAdd (R, NULL, NULL, dup,
                    W [2*p], W [2*p+1], NULL) ;
            }
            // each p writes a distinct Wnext slot from a distinct input pair
            Wnext [p] = R ;
            ownedN [p] = (R != NULL) ;

            // the two inputs are now consumed: free them immediately if we own
            // them, so only the merge frontier stays resident.  Original inputs
            // (ownedW == false) are never freed.  The unpaired trailing matrix
            // (index count-1 when count is odd) is not a pair member, so it is
            // never freed here and is carried up below.
            if (ownedW [2*p]  ) GrB_free (& W [2*p]  ) ;
            if (ownedW [2*p+1]) GrB_free (& W [2*p+1]) ;

            if (info < GrB_SUCCESS)
            {
                #pragma omp critical
                {
                    if (sum_status >= GrB_SUCCESS) sum_status = info ;
                }
            }
        }
        GRB_TRY (sum_status) ;

        // carry an unpaired (odd) trailing matrix up to the next level,
        // transferring its ownership so it is referenced (and later freed)
        // exactly once
        if (count % 2 == 1)
        {
            Wnext [npairs] = W [count-1] ;
            ownedN [npairs] = ownedW [count-1] ;
            ownedW [count-1] = false ;
        }

        // advance to the next level: swap the level buffers and their flags
        GrB_Matrix *tmpW = W ; W = Wnext ; Wnext = tmpW ;
        bool *tmpO = ownedW ; ownedW = ownedN ; ownedN = tmpO ;
        count = (count + 1) / 2 ;
    }

    //--------------------------------------------------------------------------
    // hand the root of the tree to the caller
    //--------------------------------------------------------------------------

    // With nmatrices >= 2 the final root is always a genuine sum (owned by us):
    // count==1 is only ever reached from a count==2 level, which has no carry.
    // Transfer ownership by clearing its flag so LG_FREE_WORK won't free it.
    (*C) = W [0] ;
    ownedW [0] = false ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
