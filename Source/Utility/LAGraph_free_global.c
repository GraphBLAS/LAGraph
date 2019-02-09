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

    GrB_free (&LAGraph_EQ_Complex) ;
    GrB_free (&LAGraph_SKEW_INT8) ;
    GrB_free (&LAGraph_SKEW_INT16) ;
    GrB_free (&LAGraph_SKEW_INT32) ;
    GrB_free (&LAGraph_SKEW_INT64) ;
    GrB_free (&LAGraph_SKEW_FP32) ;
    GrB_free (&LAGraph_SKEW_FP64) ;
    GrB_free (&LAGraph_SKEW_Complex) ;
    GrB_free (&LAGraph_Hermitian) ;

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

    GrB_free (&LAGraph_TRUE_BOOL) ;
    GrB_free (&LAGraph_TRUE_BOOL_Complex) ;

    GrB_free (&LAGraph_LAND_MONOID) ;
    GrB_free (&LAGraph_LOR_MONOID) ;
    GrB_free (&LAGraph_LOR_LAND_BOOL) ;

    return (GrB_SUCCESS) ;
}

