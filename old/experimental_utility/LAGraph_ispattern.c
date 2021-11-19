//------------------------------------------------------------------------------
// LAGraph_ispattern: check if a matrix is all 1
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// LAGraph_ispattern: check if a matrix values are all equal to 1.
// Contributed by Tim Davis, Texas A&M.

// FIXME: use GrB_mxv with GrB_EQ_* and bind1st

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
    GrB_free (&C) ;         \
    GrB_free (&op) ;

//****************************************************************************
// "is one" operators
#define F_UNARY(f)  ((void (*)(void *, const void *)) f)


void LAGraph_isone_int8
(
    bool *z,
    const int8_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_int16
(
    bool *z,
    const int16_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_int32
(
    bool *z,
    const int32_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_int64
(
    bool *z,
    const int64_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_uint8
(
    bool *z,
    const uint8_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_uint16
(
    bool *z,
    const uint16_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_uint32
(
    bool *z,
    const uint32_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_uint64
(
    bool *z,
    const uint64_t *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_float
(
    bool *z,
    const float *x
)
{
    (*z) = ((*x) == 1) ;
}

void LAGraph_isone_double
(
    bool *z,
    const double *x
)
{
    (*z) = ((*x) == 1) ;
}

void create_isone_op(GrB_UnaryOp *op, GrB_Type type)
{
    *op = NULL;

    if      (type == GrB_INT8  )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_int8),
                         GrB_BOOL, GrB_INT8) ;
    }
    else if (type == GrB_INT16 )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_int16),
                         GrB_BOOL, GrB_INT16) ;
    }
    else if (type == GrB_INT32 )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_int32),
                         GrB_BOOL, GrB_INT32) ;
    }
    else if (type == GrB_INT64 )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_int64),
                         GrB_BOOL, GrB_INT64) ;
    }
    else if (type == GrB_UINT8 )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_uint8),
                         GrB_BOOL, GrB_UINT8) ;
    }
    else if (type == GrB_UINT16)
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_uint16),
                         GrB_BOOL, GrB_UINT16) ;
    }
    else if (type == GrB_UINT32)
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_uint32),
                         GrB_BOOL, GrB_UINT32) ;
    }
    else if (type == GrB_UINT64)
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_uint64),
                         GrB_BOOL, GrB_UINT64) ;
    }
    else if (type == GrB_FP32  )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_float),
                         GrB_BOOL, GrB_FP32) ;
    }
    else if (type == GrB_FP64  )
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_double),
                         GrB_BOOL, GrB_FP64) ;
    }
#if 0
    else if (type == LAGraph_ComplexFP64)
    {
        GrB_UnaryOp_new (op, F_UNARY (LAGraph_isone_fc64),
                         GrB_BOOL, GxB_FC64) ;
    }
#endif
}

//****************************************************************************
GrB_Info LAGraph_ispattern  // return GrB_SUCCESS if successful
(
    bool *result,           // true if A is all one, false otherwise
    GrB_Matrix A,
    GrB_UnaryOp userop      // for A with arbitrary user-defined type.
                            // Ignored if A and B are of built-in types
)
{
#if !(LG_SUITESPARSE)
    return (GrB_PANIC) ;
#else
    GrB_Info info ;
    GrB_Matrix C = NULL ;
    GrB_UnaryOp op = NULL ;
    GrB_Type type ;
    GrB_Index nrows, ncols ;

    // check inputs
    if (result == NULL)
    {
        // error: required parameter, result, is NULL
        return (GrB_NULL_POINTER) ;
    }
    (*result) = false ;

    // get the type and size of A
    LAGRAPH_OK (GxB_Matrix_type  (&type,  A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;

    if (type == GrB_BOOL)
    {
        // result = and (A)
        LAGRAPH_OK (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, A, NULL)) ;
    }
    else
    {
        // select the unary operator
        create_isone_op(&op, type);

        if (op == NULL) op = userop ;

        if (op == NULL)
        {
            // printf ("LAGraph_ispattern: userop is NULL\n") ;
            return (GrB_NULL_POINTER) ;
        }

        // C = isone (A)
        LAGRAPH_OK (GrB_Matrix_new (&C, GrB_BOOL, nrows, ncols)) ;

        LAGRAPH_OK (GrB_apply (C, NULL, NULL, op, A, NULL)) ;

        // result = and (C)
        LAGRAPH_OK (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, C, NULL)) ;
    }

    // free workspace and return result
    LAGraph_FREE_ALL ;
    return (GrB_SUCCESS) ;
#endif
}
