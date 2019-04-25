//------------------------------------------------------------------------------
// LAGraph_free_global:  free all global operators and types
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Free all global operators and types for LAGraph.

#include "LAGraph_internal.h"

GrB_Info LAGraph_free_global ( )
{

    // free the types
    GrB_free (&LAGraph_Complex) ;

    // free the unary ops
    GrB_free (&LAGraph_ISONE_INT8) ;
    GrB_free (&LAGraph_ISONE_INT16) ;
    GrB_free (&LAGraph_ISONE_INT32) ;
    GrB_free (&LAGraph_ISONE_INT64) ;
    GrB_free (&LAGraph_ISONE_UINT8) ;
    GrB_free (&LAGraph_ISONE_UINT16) ;
    GrB_free (&LAGraph_ISONE_UINT32) ;
    GrB_free (&LAGraph_ISONE_UINT64) ;
    GrB_free (&LAGraph_ISONE_FP32) ;
    GrB_free (&LAGraph_ISONE_FP64) ;
    GrB_free (&LAGraph_ISONE_Complex) ;
    GrB_free (&LAGraph_ISTWO_UINT32) ;
    GrB_free (&LAGraph_DECR_INT32) ;
    GrB_free (&LAGraph_DECR_INT64) ;
    GrB_free (&LAGraph_TRUE_BOOL) ;
    GrB_free (&LAGraph_TRUE_BOOL_Complex) ;
    GrB_free (&LAGraph_ONE_FP64) ;
    GrB_free (&LAGraph_ONE_UINT32) ;
    GrB_free (&LAGraph_COMB_FP64) ;

    // free the binary ops
    GrB_free (&LAGraph_EQ_Complex) ;
    GrB_free (&LAGraph_SKEW_INT8) ;
    GrB_free (&LAGraph_SKEW_INT16) ;
    GrB_free (&LAGraph_SKEW_INT32) ;
    GrB_free (&LAGraph_SKEW_INT64) ;
    GrB_free (&LAGraph_SKEW_FP32) ;
    GrB_free (&LAGraph_SKEW_FP64) ;
    GrB_free (&LAGraph_SKEW_Complex) ;
    GrB_free (&LAGraph_Hermitian) ;
    GrB_free (&LAGraph_LOR_UINT32) ;

    // free the monoids
    GrB_free (&LAGraph_PLUS_INT64_MONOID) ;
    GrB_free (&LAGraph_MAX_INT32_MONOID) ;
    GrB_free (&LAGraph_LAND_MONOID) ;
    GrB_free (&LAGraph_LOR_MONOID) ;
    GrB_free (&LAGraph_MIN_INT32_MONOID) ;
    GrB_free (&LAGraph_MIN_INT64_MONOID) ;
    GrB_free (&LAGraph_PLUS_UINT32_MONOID) ;
    GrB_free (&LAGraph_PLUS_FP64_MONOID) ;

    // free the semirings
    GrB_free (&LAGraph_LOR_LAND_BOOL) ;
    GrB_free (&LAGraph_LOR_SECOND_BOOL) ;
    GrB_free (&LAGraph_LOR_FIRST_BOOL) ;
    GrB_free (&LAGraph_MIN_SECOND_INT32) ;
    GrB_free (&LAGraph_MIN_FIRST_INT32) ;
    GrB_free (&LAGraph_MIN_SECOND_INT64) ;
    GrB_free (&LAGraph_MIN_FIRST_INT64) ;
    GrB_free (&LAGraph_PLUS_TIMES_UINT32) ;
    GrB_free (&LAGraph_PLUS_TIMES_FP64) ;

    // free the descriptors
//  GrB_free (&LAGraph_desc_oooo) ;     // NULL, no need to free it
    GrB_free (&LAGraph_desc_ooor) ;
    GrB_free (&LAGraph_desc_ooco) ;
    GrB_free (&LAGraph_desc_oocr) ;
    GrB_free (&LAGraph_desc_otoo) ;
    GrB_free (&LAGraph_desc_otor) ;
    GrB_free (&LAGraph_desc_otco) ;
    GrB_free (&LAGraph_desc_otcr) ;
    GrB_free (&LAGraph_desc_tooo) ;
    GrB_free (&LAGraph_desc_toor) ;
    GrB_free (&LAGraph_desc_toco) ;
    GrB_free (&LAGraph_desc_tocr) ;
    GrB_free (&LAGraph_desc_ttoo) ;
    GrB_free (&LAGraph_desc_ttor) ;
    GrB_free (&LAGraph_desc_ttco) ;
    GrB_free (&LAGraph_desc_ttcr) ;

    // free the select operator for ktruss and allktruss
    GrB_free (&LAGraph_support) ;

    return (GrB_SUCCESS) ;
}

