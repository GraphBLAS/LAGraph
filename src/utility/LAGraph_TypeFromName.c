//------------------------------------------------------------------------------
// LAGraph_TypeFromName: return the GrB_Type corresponding to its given name
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2023 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// This method works for any GraphBLAS library.  On input, name is a char array
// of length at least LAGRAPH_MAX_NAME_LEN.

// Only built-in types are supported.  User-defined types are not supported.

#include "LG_internal.h"

int LAGraph_TypeFromName
(
    // output:
    GrB_Type *type, // GraphBLAS type
    // input:
    char *name,     // name of the type: a null-terminated string
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (type != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the GrB_Type from its name
    //--------------------------------------------------------------------------

    #define MATCH2(s1,s2) MATCHNAME (name, s1) || MATCHNAME (name, s2)

    if      (MATCH2 ("bool"    , "GrB_BOOL"   )) (*type) = GrB_BOOL   ;
    else if (MATCH2 ("int8_t"  , "GrB_INT8"   )) (*type) = GrB_INT8   ;
    else if (MATCH2 ("int16_t" , "GrB_INT16"  )) (*type) = GrB_INT16  ;
    else if (MATCH2 ("int32_t" , "GrB_INT32"  )) (*type) = GrB_INT32  ;
    else if (MATCH2 ("int64_t" , "GrB_INT64"  )) (*type) = GrB_INT64  ;
    else if (MATCH2 ("uint8_t" , "GrB_UINT8"  )) (*type) = GrB_UINT8  ;
    else if (MATCH2 ("uint16_t", "GrB_UINT16" )) (*type) = GrB_UINT16 ;
    else if (MATCH2 ("uint32_t", "GrB_UINT32" )) (*type) = GrB_UINT32 ;
    else if (MATCH2 ("uint64_t", "GrB_UINT64" )) (*type) = GrB_UINT64 ;
    else if (MATCH2 ("float"   , "GrB_FP32"   )) (*type) = GrB_FP32   ;
    else if (MATCH2 ("double"  , "GrB_FP64"   )) (*type) = GrB_FP64   ;
    // if complex types from SuiteSparse:GraphBLAS are added to LAGraph:
//  else if (MATCH2 ("float  complex", "GxB_FC32") ||
//           MATCHNAME (name, "float _Complex" )) (*type) = GxB_FC32   ;
//  else if (MATCH2 ("double complex", "GxB_FC64") ||
//           MATCHNAME (name, "double _Complex")) (*type) = GxB_FC64   ;
    else
    {
        (*type) = NULL ;
        LG_ASSERT_MSGF (false, GrB_NOT_IMPLEMENTED,
            "type \"%s\" not supported", name) ;
    }
    return (GrB_SUCCESS) ;
}

