//------------------------------------------------------------------------------
// LAGraph_TypeName: return the name of a type
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_Internal.h"

int LAGraph_TypeName        // returns 0 if successful, -1 if failure
(
    char **name,            // name of the type
    GrB_Type type,          // GraphBLAS type
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LAGraph_CLEAR_MSG ;
    LAGraph_CHECK (type == NULL, -1, "type is NULL") ;
    LAGraph_CHECK (name == NULL, -1, "name is NULL") ;

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
    else if (type == GxB_FC32  ) (*name) = "single complex" ;
    else if (type == GxB_FC64  ) (*name) = "double complex" ;
    else (*name) = "user-defined" ;

    return (0) ;
}

