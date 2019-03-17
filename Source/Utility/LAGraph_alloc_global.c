//------------------------------------------------------------------------------
// LAGraph_alloc_global:  allocate all global objects for LAGraph
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
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
    // printf ("eq complex:\nx (%.18e,%.18e)\ny (%.18e,%.18e)\n",
    //     creal (*x), cimag (*x), creal (*y), cimag (*y)) ;
    (*z) = (*x) == (*y) ;
    // printf ("result %d\n", (*z)) ;
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

// unary operators that return boolean true
GrB_UnaryOp LAGraph_TRUE_BOOL = NULL ;
GrB_UnaryOp LAGraph_TRUE_BOOL_Complex  = NULL ;

void LAGraph_true_bool
(
    bool *z,
    const bool *x       // ignored
)
{
    (*z) = true ;
}

void LAGraph_true_bool_complex
(
    bool *z,
    const double complex *x     // ignored
)
{
    (*z) = true ;
}

// integer decrement
GrB_UnaryOp LAGraph_DECR_INT32 = NULL ;
GrB_UnaryOp LAGraph_DECR_INT64 = NULL ;

void LAGraph_decr_int32
(
    int32_t *z,
    const int32_t *x
)
{
    (*z) = ((*x) - 1) ;
}

void LAGraph_decr_int64
(
    int64_t *z,
    const int64_t *x
)
{
    (*z) = ((*x) - 1) ;
}

// boolean monoids and semirings
GrB_Monoid LAGraph_PLUS_INT64_MONOID = NULL ;
GrB_Monoid LAGraph_MAX_INT32_MONOID = NULL ;
GrB_Monoid LAGraph_LAND_MONOID = NULL ;
GrB_Monoid LAGraph_LOR_MONOID = NULL ;

GrB_Monoid LAGraph_MIN_INT32_MONOID = NULL ;
GrB_Monoid LAGraph_MIN_INT64_MONOID = NULL ;

GrB_Semiring LAGraph_LOR_LAND_BOOL = NULL ;

GrB_Semiring LAGraph_LOR_SECOND_BOOL = NULL ;
GrB_Semiring LAGraph_LOR_FIRST_BOOL = NULL ;

GrB_Semiring LAGraph_MIN_SECOND_INT32 = NULL ;
GrB_Semiring LAGraph_MIN_FIRST_INT32 = NULL ;

GrB_Semiring LAGraph_MIN_SECOND_INT64 = NULL ;
GrB_Semiring LAGraph_MIN_FIRST_INT64 = NULL ;

// all 16 descriptors
// syntax: 4 characters define the following.  'o' is the default:
// 1: o or t: A transpose
// 2: o or t: B transpose
// 3: o or c: complemented mask
// 4: o or r: replace
GrB_Descriptor

    LAGraph_desc_oooo = NULL ,   // default (NULL)
    LAGraph_desc_ooor = NULL ,   // replace
    LAGraph_desc_ooco = NULL ,   // compl mask
    LAGraph_desc_oocr = NULL ,   // compl mask, replace

    LAGraph_desc_otoo = NULL ,   // A'
    LAGraph_desc_otor = NULL ,   // A', replace
    LAGraph_desc_otco = NULL ,   // A', compl mask
    LAGraph_desc_otcr = NULL ,   // A', compl mask, replace

    LAGraph_desc_tooo = NULL ,   // B'
    LAGraph_desc_toor = NULL ,   // B', replace
    LAGraph_desc_toco = NULL ,   // B', compl mask
    LAGraph_desc_tocr = NULL ,   // B', compl mask, replace

    LAGraph_desc_ttoo = NULL ,   // A', B'
    LAGraph_desc_ttor = NULL ,   // A', B', replace
    LAGraph_desc_ttco = NULL ,   // A', B', compl mask
    LAGraph_desc_ttcr = NULL ;   // A', B', compl mask, replace

//------------------------------------------------------------------------------
// LAGraph_alloc_global
//------------------------------------------------------------------------------

#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)
#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

GrB_Info LAGraph_alloc_global ( )
{
    GrB_Info info ;

    //--------------------------------------------------------------------------
    // create the complex type for LAGraph
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Type_new (&LAGraph_Complex, sizeof (double complex))) ;

    //--------------------------------------------------------------------------
    // create the binary operators
    //--------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------
    // create the unary operators that check if equal to 1
    //--------------------------------------------------------------------------

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

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_FP32,
        F_UNARY (LAGraph_isone_float),
        GrB_BOOL, GrB_FP32)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_FP64,
        F_UNARY (LAGraph_isone_double),
        GrB_BOOL, GrB_FP64)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ISONE_Complex,
        F_UNARY (LAGraph_isone_complex),
        GrB_BOOL, LAGraph_Complex)) ;

    //--------------------------------------------------------------------------
    // create the unary decrement operators
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_DECR_INT32,
        F_UNARY (LAGraph_decr_int32),
        GrB_INT32, GrB_INT32)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_DECR_INT64,
        F_UNARY (LAGraph_decr_int64),
        GrB_INT64, GrB_INT64)) ;

    //--------------------------------------------------------------------------
    // create the unary operators that return true
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_TRUE_BOOL,
        F_UNARY (LAGraph_true_bool),
        GrB_BOOL, GrB_BOOL)) ;

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_TRUE_BOOL_Complex,
        F_UNARY (LAGraph_true_bool_complex),
        GrB_BOOL, LAGraph_Complex)) ;

    //--------------------------------------------------------------------------
    // create the monoids
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Monoid_new_INT64 (&LAGraph_PLUS_INT64_MONOID,
        GrB_PLUS_INT64, 0)) ;

    LAGRAPH_OK (GrB_Monoid_new_INT32 (&LAGraph_MAX_INT32_MONOID,
        GrB_MAX_INT32, INT32_MIN)) ;

    LAGRAPH_OK (GrB_Monoid_new_INT32 (&LAGraph_MIN_INT32_MONOID,
        GrB_MIN_INT32, INT32_MAX)) ;

    LAGRAPH_OK (GrB_Monoid_new_INT64 (&LAGraph_MIN_INT64_MONOID,
        GrB_MIN_INT64, INT64_MAX)) ;

    LAGRAPH_OK (GrB_Monoid_new_BOOL (&LAGraph_LAND_MONOID, GrB_LAND, true )) ;
    LAGRAPH_OK (GrB_Monoid_new_BOOL (&LAGraph_LOR_MONOID , GrB_LOR , false)) ;

    //--------------------------------------------------------------------------
    // create the semirings
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_LOR_LAND_BOOL,
        LAGraph_LOR_MONOID, GrB_LAND)) ;

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_LOR_FIRST_BOOL,
        LAGraph_LOR_MONOID, GrB_FIRST_BOOL)) ;

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_LOR_SECOND_BOOL,
        LAGraph_LOR_MONOID, GrB_SECOND_BOOL)) ;

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_MIN_SECOND_INT32,
        LAGraph_MIN_INT32_MONOID, GrB_SECOND_INT32)) ;

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_MIN_FIRST_INT32,
        LAGraph_MIN_INT32_MONOID, GrB_FIRST_INT32)) ;

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_MIN_SECOND_INT64,
        LAGraph_MIN_INT64_MONOID, GrB_SECOND_INT64)) ;

    LAGRAPH_OK (GrB_Semiring_new (&LAGraph_MIN_FIRST_INT64,
        LAGraph_MIN_INT64_MONOID, GrB_FIRST_INT64)) ;

    //--------------------------------------------------------------------------
    // create 15 descriptors (one does not need to be allocated)
    //--------------------------------------------------------------------------

    LAGraph_desc_oooo = NULL ;   // default (NULL)
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_ooor)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_ooco)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_oocr)) ;

    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_otoo)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_otor)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_otco)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_otcr)) ;

    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_tooo)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_toor)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_toco)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_tocr)) ;

    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_ttoo)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_ttor)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_ttco)) ;
    LAGRAPH_OK (GrB_Descriptor_new (&LAGraph_desc_ttcr)) ;

    //--------------------------------------------------------------------------
    // set the descriptors
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ooor, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_oocr, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otor, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otcr, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_toor, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_tocr, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttor, GrB_OUTP, GrB_REPLACE)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttcr, GrB_OUTP, GrB_REPLACE)) ;

    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ooco, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_oocr, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otco, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otcr, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_toco, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_tocr, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttco, GrB_MASK, GrB_SCMP)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttcr, GrB_MASK, GrB_SCMP)) ;

    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otoo, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otor, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otco, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_otcr, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttoo, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttor, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttco, GrB_INP1, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttcr, GrB_INP1, GrB_TRAN)) ;

    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_tooo, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_toor, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_toco, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_tocr, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttoo, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttor, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttco, GrB_INP0, GrB_TRAN)) ;
    LAGRAPH_OK (GrB_Descriptor_set (LAGraph_desc_ttcr, GrB_INP0, GrB_TRAN)) ;

    return (GrB_SUCCESS) ;
}

