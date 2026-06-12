//------------------------------------------------------------------------------
// LAGraph_Matrix_Sum: sum an array of matrices with a binary operator
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

// LAGraph_Matrix_Sum combines an array of matrices into a single matrix C.  It
// computes the total number of entries across all inputs, allocates a single
// tuple buffer (I, J, X) large enough to hold every entry, extracts the tuples
// of each input matrix into that buffer, and then calls GrB_Matrix_build with
// the binary operator dup to combine any duplicate (i,j) entries.  With dup =
// GrB_PLUS_FP64 (for example) this computes the element-wise sum of all input
// matrices.

// All input matrices must have identical dimensions and identical built-in
// type; C is created with that same type and dimensions.

#define LG_FREE_WORK                            \
{                                               \
    LAGraph_Free ((void **) &I, NULL) ;         \
    LAGraph_Free ((void **) &J, NULL) ;         \
    LAGraph_Free ((void **) &X, NULL) ;         \
}

#define LG_FREE_ALL                             \
{                                               \
    LG_FREE_WORK ;                              \
    GrB_free (C) ;                              \
}

#include "LG_internal.h"

int LAGraph_Matrix_Sum
(
    // output:
    GrB_Matrix *C,          // result = combination of all input matrices
    // input:
    GrB_Matrix *Matrices,   // array of nmatrices input matrices
    GrB_Index nmatrices,    // number of matrices in the array (must be >= 1)
    GrB_BinaryOp dup,       // operator to combine duplicate (i,j) entries
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Index *I = NULL, *J = NULL ;
    void *X = NULL ;
    LG_ASSERT_MSG (C != NULL, GrB_NULL_POINTER, "&C != NULL") ;
    LG_ASSERT (Matrices != NULL, GrB_NULL_POINTER) ;
    (*C) = NULL ;
    LG_ASSERT_MSG (nmatrices >= 1, GrB_INVALID_VALUE,
        "nmatrices must be >= 1") ;
    LG_ASSERT (Matrices [0] != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the reference dimensions and type from the first matrix
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    int32_t typecode ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, Matrices [0])) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, Matrices [0])) ;
    GRB_TRY (GrB_get (Matrices [0], &typecode, GrB_EL_TYPE_CODE)) ;

    //--------------------------------------------------------------------------
    // validate every matrix and accumulate the total number of entries
    //--------------------------------------------------------------------------

    GrB_Index total = 0 ;
    for (GrB_Index k = 0 ; k < nmatrices ; k++)
    {
        GrB_Matrix Ak = Matrices [k] ;
        LG_ASSERT (Ak != NULL, GrB_NULL_POINTER) ;
        GrB_Index r, c, n ;
        int32_t code ;
        GRB_TRY (GrB_Matrix_nrows (&r, Ak)) ;
        GRB_TRY (GrB_Matrix_ncols (&c, Ak)) ;
        LG_ASSERT_MSG (r == nrows && c == ncols, GrB_DIMENSION_MISMATCH,
            "all input matrices must have the same dimensions") ;
        GRB_TRY (GrB_get (Ak, &code, GrB_EL_TYPE_CODE)) ;
        LG_ASSERT_MSG (code == typecode, GrB_DOMAIN_MISMATCH,
            "all input matrices must have the same type") ;
        GRB_TRY (GrB_Matrix_nvals (&n, Ak)) ;
        total += n ;
    }

    //--------------------------------------------------------------------------
    // allocate the shared row/column index buffers (guard against size 0)
    //--------------------------------------------------------------------------

    GrB_Index alloc = (total == 0) ? 1 : total ;
    LG_TRY (LAGraph_Malloc ((void **) &I, alloc, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &J, alloc, sizeof (GrB_Index), msg)) ;

    //--------------------------------------------------------------------------
    // extract tuples from every matrix, then build the result
    //--------------------------------------------------------------------------

    // For each built-in type: allocate the value buffer X with the correct
    // element size, extract the tuples of every input matrix into the shared
    // buffer at the running offset, create C, and build it with the dup
    // operator to combine duplicate (i,j) entries.

    #define LG_SUM_CASE(code, ctype, gtype, suffix)                          \
        case code :                                                          \
        {                                                                    \
            ctype *Xt = NULL ;                                               \
            LG_TRY (LAGraph_Malloc ((void **) &Xt, alloc, sizeof (ctype),    \
                msg)) ;                                                      \
            X = (void *) Xt ;                                                \
            GrB_Index offset = 0 ;                                           \
            for (GrB_Index k = 0 ; k < nmatrices ; k++)                      \
            {                                                                \
                GrB_Index n, got ;                                           \
                GRB_TRY (GrB_Matrix_nvals (&n, Matrices [k])) ;              \
                if (n == 0) continue ;                                       \
                got = n ;                                                    \
                GRB_TRY (GrB_Matrix_extractTuples_ ## suffix (               \
                    I + offset, J + offset, Xt + offset, &got,               \
                    Matrices [k])) ;                                         \
                offset += n ;                                                \
            }                                                                \
            GRB_TRY (GrB_Matrix_new (C, gtype, nrows, ncols)) ;              \
            GRB_TRY (GrB_Matrix_build_ ## suffix (*C, I, J, Xt, total,       \
                dup)) ;                                                      \
        }                                                                    \
        break ;

    switch (typecode)
    {
        LG_SUM_CASE (GrB_BOOL_CODE,   bool,     GrB_BOOL,   BOOL  )
        LG_SUM_CASE (GrB_INT8_CODE,   int8_t,   GrB_INT8,   INT8  )
        LG_SUM_CASE (GrB_INT16_CODE,  int16_t,  GrB_INT16,  INT16 )
        LG_SUM_CASE (GrB_INT32_CODE,  int32_t,  GrB_INT32,  INT32 )
        LG_SUM_CASE (GrB_INT64_CODE,  int64_t,  GrB_INT64,  INT64 )
        LG_SUM_CASE (GrB_UINT8_CODE,  uint8_t,  GrB_UINT8,  UINT8 )
        LG_SUM_CASE (GrB_UINT16_CODE, uint16_t, GrB_UINT16, UINT16)
        LG_SUM_CASE (GrB_UINT32_CODE, uint32_t, GrB_UINT32, UINT32)
        LG_SUM_CASE (GrB_UINT64_CODE, uint64_t, GrB_UINT64, UINT64)
        LG_SUM_CASE (GrB_FP32_CODE,   float,    GrB_FP32,   FP32  )
        LG_SUM_CASE (GrB_FP64_CODE,   double,   GrB_FP64,   FP64  )
        default :
            LG_ASSERT_MSG (false, GrB_NOT_IMPLEMENTED,
                "user-defined types are not supported") ;
    }

    #undef LG_SUM_CASE

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
