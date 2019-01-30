//------------------------------------------------------------------------------
// LAGraph_alloc_global:  allocate all global objects for LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Define and allocate global types and operators for LAGraph.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL                    \
{                                           \
    LAGraph_free_global ( ) ;               \
}

// a global value for returning the complex type in a Matrix Market file:
GrB_Type LAGraph_Complex = NULL ;

// binary operators to test for symmetry, skew-symmetry and Hermitian property
GrB_BinaryOp LAGraph_EQ_Complex   = NULL ;
GrB_BinaryOp LAGraph_SKEW_INT8    = NULL ;
GrB_BinaryOp LAGraph_SKEW_INT16   = NULL ;
GrB_BinaryOp LAGraph_SKEW_INT32   = NULL ;
GrB_BinaryOp LAGraph_SKEW_INT64   = NULL ;
GrB_BinaryOp LAGraph_SKEW_FP32    = NULL ;
GrB_BinaryOp LAGraph_SKEW_FP64    = NULL ;
GrB_BinaryOp LAGraph_SKEW_Complex = NULL ;
GrB_BinaryOp LAGraph_Hermitian    = NULL ;

void LAGraph_eq_complex
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = (*x) == (*y) ;
}

void LAGraph_skew_int8
(
    bool *z,
    const int8_t *x,
    const int8_t *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_skew_int16
(
    bool *z,
    const int16_t *x,
    const int16_t *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_skew_int32
(
    bool *z,
    const int32_t *x,
    const int32_t *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_skew_int64
(
    bool *z,
    const int64_t *x,
    const int64_t *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_skew_float
(
    bool *z,
    const float *x,
    const float *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_skew_double
(
    bool *z,
    const double *x,
    const double *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_skew_complex
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = (*x) == -(*y) ;
}

void LAGraph_hermitian
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    (*z) = (*x) == conj (*y) ;
}

// unary operators to check if the entry is equal to 1
GrB_UnaryOp LAGraph_ISONE_INT8     = NULL ;
GrB_UnaryOp LAGraph_ISONE_INT16    = NULL ;
GrB_UnaryOp LAGraph_ISONE_INT32    = NULL ;
GrB_UnaryOp LAGraph_ISONE_INT64    = NULL ;
GrB_UnaryOp LAGraph_ISONE_UINT8    = NULL ;
GrB_UnaryOp LAGraph_ISONE_UINT16   = NULL ;
GrB_UnaryOp LAGraph_ISONE_UINT32   = NULL ;
GrB_UnaryOp LAGraph_ISONE_UINT64   = NULL ;
GrB_UnaryOp LAGraph_ISONE_FP32     = NULL ;
GrB_UnaryOp LAGraph_ISONE_FP64     = NULL ;
GrB_UnaryOp LAGraph_ISONE_Complex  = NULL ;

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

void LAGraph_isone_complex
(
    bool *z,
    const double complex *x
)
{
    (*z) = ((*x) == 1) ;
}

// boolean monoid
GrB_Monoid LAGraph_LAND_MONOID = NULL ;

//------------------------------------------------------------------------------
// LAGraph_alloc_global
//------------------------------------------------------------------------------

#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)
#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

GrB_Info LAGraph_alloc_global ( )
{

    // create the complex type for LAGraph

    LAGRAPH_OK (GrB_Type_new (&LAGraph_Complex, sizeof (double complex))) ;

    // create the binary operators

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_EQ_Complex,
        F_BINARY (LAGraph_eq_complex),
        GrB_BOOL, LAGraph_Complex, LAGraph_Complex)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_INT8,
        F_BINARY (LAGraph_skew_int8),
        GrB_BOOL, GrB_INT8, GrB_INT8)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_INT16,
        F_BINARY (LAGraph_skew_int16),
        GrB_BOOL, GrB_INT16, GrB_INT16)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_INT32,
        F_BINARY (LAGraph_skew_int32),
        GrB_BOOL, GrB_INT32, GrB_INT32)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_INT64, 
        F_BINARY (LAGraph_skew_int64),
        GrB_BOOL, GrB_INT64, GrB_INT64)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_FP32, 
        F_BINARY (LAGraph_skew_float),
        GrB_BOOL, GrB_FP32, GrB_FP32)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_FP64, 
        F_BINARY (LAGraph_skew_double),
        GrB_BOOL, GrB_FP64, GrB_FP64)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_SKEW_Complex, 
        F_BINARY (LAGraph_skew_complex),
        GrB_BOOL, LAGraph_Complex, LAGraph_Complex)) ;

    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_Hermitian, 
        F_BINARY (LAGraph_hermitian),
        GrB_BOOL, LAGraph_Complex, LAGraph_Complex)) ;

    // create the unary operators

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_INT8, 
        F_UNARY (LAGraph_isone_int8),
        GrB_BOOL, GrB_INT8)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_INT16, 
        F_UNARY (LAGraph_isone_int16),
        GrB_BOOL, GrB_INT16)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_INT32, 
        F_UNARY (LAGraph_isone_int32),
        GrB_BOOL, GrB_INT32)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_INT64, 
        F_UNARY (LAGraph_isone_int64),
        GrB_BOOL, GrB_INT64)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT8, 
        F_UNARY (LAGraph_isone_uint8),
        GrB_BOOL, GrB_UINT8)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT16, 
        F_UNARY (LAGraph_isone_uint16),
        GrB_BOOL, GrB_UINT16)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT32, 
        F_UNARY (LAGraph_isone_uint32),
        GrB_BOOL, GrB_UINT32)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT64, 
        F_UNARY (LAGraph_isone_uint64),
        GrB_BOOL, GrB_UINT64)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT64, 
        F_UNARY (LAGraph_isone_float),
        GrB_BOOL, GrB_FP32)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT64, 
        F_UNARY (LAGraph_isone_double),
        GrB_BOOL, GrB_FP64)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_UINT64, 
        F_UNARY (LAGraph_isone_complex),
        GrB_BOOL, LAGraph_Complex)) ;

    LAGRAPH_OK (GrB_Monoid_new_BOOL (&LAGraph_LAND_MONOID, GrB_LAND, true)) ;

    return (GrB_SUCCESS) ;
}

