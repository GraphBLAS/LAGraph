//------------------------------------------------------------------------------
// LAGraph_TypeName: return the name of a type
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_TypeName        // returns 0 if successful, < 0 if failure
(
    char **name,            // name of the type
    GrB_Type type,          // GraphBLAS type
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (type == NULL, GrB_NULL_POINTER, "type is NULL") ;
    LG_CHECK (name == NULL, GrB_NULL_POINTER, "name is NULL") ;

    //--------------------------------------------------------------------------
    // determine the name of the type
    //--------------------------------------------------------------------------

    if      (type == GrB_BOOL  ) (*name) = "bool"   ;
    else if (type == GrB_INT8  ) (*name) = "int8"   ;
    else if (type == GrB_INT16 ) (*name) = "int16"  ;
    else if (type == GrB_INT32 ) (*name) = "int32"  ;
    else if (type == GrB_INT64 ) (*name) = "int64"  ;
    else if (type == GrB_UINT8 ) (*name) = "uint8"  ;
    else if (type == GrB_UINT16) (*name) = "uint16" ;
    else if (type == GrB_UINT32) (*name) = "uint32" ;
    else if (type == GrB_UINT64) (*name) = "uint64" ;
    else if (type == GrB_FP32  ) (*name) = "single" ;
    else if (type == GrB_FP64  ) (*name) = "double" ;
    #if 0
    else if (type == GxB_FC32  ) (*name) = "single complex" ;
    else if (type == GxB_FC64  ) (*name) = "double complex" ;
    #endif
    else (*name) = "user-defined" ;

    return (0) ;
}

