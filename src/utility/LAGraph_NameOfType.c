//------------------------------------------------------------------------------
// LAGraph_NameOfType: return the C name of a GraphBLAS GrB_Type
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

#define LG_FREE_WORK LAGraph_Free ((void **) &s, NULL) ;
#define LG_FREE_ALL LG_FREE_WORK

#include "LG_internal.h"

int LAGraph_NameOfType
(
    // output:
    char *name,     // name of the type: user provided array of size at
                    // least LAGRAPH_MAX_NAME_LEN.
    // input:
    GrB_Type type,  // GraphBLAS type
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    char *s = NULL ;
    LG_CLEAR_MSG ;
    LG_ASSERT (type != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the name of the type
    //--------------------------------------------------------------------------

    size_t len ;
    name [0] = '\0' ;

    int32_t typecode ;
    GRB_TRY (GrB_Type_get_INT32 (type, &typecode, GrB_EL_TYPE_CODE)) ;

    switch (typecode)
    {

        // for user-defined types, return the GrB_NAME of the type
        default :
        case GrB_UDT_CODE    : 
            GRB_TRY (GrB_Type_get_SIZE (type, &len, GrB_NAME)) ;
            LG_TRY (LAGraph_Malloc ((void **) &s, len+1, sizeof (char), msg)) ;
            GRB_TRY (GrB_Type_get_String (type, s, GrB_NAME)) ;
            len = LAGRAPH_MIN (len, LAGRAPH_MAX_NAME_LEN) ;
            strncpy (name, s, len) ;
            name [LAGRAPH_MAX_NAME_LEN-1] = '\0' ;
            LG_FREE_WORK ;
            break ;

        // for built-in types, return the C type name
        case GrB_BOOL_CODE   : strcpy (name, "bool"    ) ; break ;
        case GrB_INT8_CODE   : strcpy (name, "int8_t"  ) ; break ;
        case GrB_INT16_CODE  : strcpy (name, "int16_t" ) ; break ;
        case GrB_INT32_CODE  : strcpy (name, "int32_t" ) ; break ;
        case GrB_INT64_CODE  : strcpy (name, "int64_t" ) ; break ;
        case GrB_UINT8_CODE  : strcpy (name, "uint8_t" ) ; break ;
        case GrB_UINT16_CODE : strcpy (name, "uint16_t") ; break ;
        case GrB_UINT32_CODE : strcpy (name, "uint32_t") ; break ;
        case GrB_UINT64_CODE : strcpy (name, "uint64_t") ; break ;
        case GrB_FP32_CODE   : strcpy (name, "float"   ) ; break ;
        case GrB_FP64_CODE   : strcpy (name, "double"  ) ; break ;
//      case GxB_FC32_CODE   : strcpy (name, "float complex"  ) ; break ;
//      case GxB_FC64_CODE   : strcpy (name, "double complex" ) ; break ;
    }
    return (GrB_SUCCESS) ;
}

