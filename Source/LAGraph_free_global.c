//------------------------------------------------------------------------------
// LAGraph_free_global:  free all global operators and types
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

GrB_Info LAGraph_free_global ( )
{

    GrB_free (&LAGraph_Complex) ;
    GrB_free (&LAGraph_SKEW_INT8) ;
    GrB_free (&LAGraph_SKEW_INT16) ;
    GrB_free (&LAGraph_SKEW_INT32) ;
    GrB_free (&LAGraph_SKEW_INT64) ;
    GrB_free (&LAGraph_SKEW_FP32) ;
    GrB_free (&LAGraph_SKEW_FP64) ;
    GrB_free (&LAGraph_SKEW_Complex) ;
    GrB_free (&LAGraph_Hermitian) ;

    return (GrB_SUCCESS) ;
}

