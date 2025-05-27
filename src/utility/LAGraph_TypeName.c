//------------------------------------------------------------------------------
// LAGraph_*TypeName: return the name of type of a matrix, vector, or scalar
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

// On input, "char *name" is a pointer to a pre-allocated array of size at
// least LAGRAPH_MAX_NAME_LEN.  On output, the array is filled with a string
// corresponding to the type of a GrB_Matrix, GrB_Vector, or GrB_Scalar.
// For built-in types, the strings are defined as:

//      "bool"      GrB_BOOL
//      "int8_t"    GrB_INT8
//      "int16_t"   GrB_INT16
//      "int32_t"   GrB_INT32
//      "int64_t"   GrB_INT64
//      "uint8_t"   GrB_UINT8
//      "uint16_t"  GrB_UINT16
//      "uint32_t"  GrB_UINT32
//      "uint64_t"  GrB_UINT64
//      "float"     GrB_FP32
//      "double"    GrB_FP64

// For user-defined types, the GrB_NAME of the type is returned.

#define LG_FREE_WORK LAGraph_Free ((void **) &t, NULL) ;
#define LG_FREE_ALL LG_FREE_WORK

#include "LG_internal.h"

#define TYPENAME(object)                                                       \
    char *t = NULL ;                                                           \
    LG_CLEAR_MSG ;                                                             \
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;                               \
    size_t len ;                                                               \
    name [0] = '\0' ;                                                          \
    int32_t typecode ;                                                         \
    GRB_TRY (GrB_get (object, &typecode, GrB_EL_TYPE_CODE)) ;                  \
    switch (typecode)                                                          \
    {                                                                          \
        /* for user-defined types, return the GrB_EL_TYPE_STRING */            \
        default :                                                              \
        case GrB_UDT_CODE    :                                                 \
            GRB_TRY (GrB_get (object, &len, GrB_EL_TYPE_STRING)) ;             \
            LG_TRY (LAGraph_Malloc ((void **) &t, len+1, sizeof (char), msg)) ;\
            GRB_TRY (GrB_get (object, t, GrB_EL_TYPE_STRING)) ;                \
            len = LAGRAPH_MIN (len, LAGRAPH_MAX_NAME_LEN) ;                    \
            strncpy (name, t, len) ;                                           \
            name [LAGRAPH_MAX_NAME_LEN-1] = '\0' ;                             \
            LG_FREE_WORK ;                                                     \
            break ;                                                            \
        /* for built-in types, return the C type name */                       \
        case GrB_BOOL_CODE   : strcpy (name, "bool"    ) ; break ;             \
        case GrB_INT8_CODE   : strcpy (name, "int8_t"  ) ; break ;             \
        case GrB_INT16_CODE  : strcpy (name, "int16_t" ) ; break ;             \
        case GrB_INT32_CODE  : strcpy (name, "int32_t" ) ; break ;             \
        case GrB_INT64_CODE  : strcpy (name, "int64_t" ) ; break ;             \
        case GrB_UINT8_CODE  : strcpy (name, "uint8_t" ) ; break ;             \
        case GrB_UINT16_CODE : strcpy (name, "uint16_t") ; break ;             \
        case GrB_UINT32_CODE : strcpy (name, "uint32_t") ; break ;             \
        case GrB_UINT64_CODE : strcpy (name, "uint64_t") ; break ;             \
        case GrB_FP32_CODE   : strcpy (name, "float"   ) ; break ;             \
        case GrB_FP64_CODE   : strcpy (name, "double"  ) ; break ;             \
/*      case GxB_FC32_CODE   : strcpy (name, "float complex"  ) ; break ; */   \
/*      case GxB_FC64_CODE   : strcpy (name, "double complex" ) ; break ; */   \
    }                                                                          \
    return (GrB_SUCCESS) ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_TypeName: return the name of the GrB_Type of a GrB_Matrix
//------------------------------------------------------------------------------

int LAGraph_Matrix_TypeName
(
    // output:
    char *name,     // name of the type of the matrix A (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    // input:
    GrB_Matrix A,   // matrix to query
    char *msg
)
{
    TYPENAME (A) ;
}

//------------------------------------------------------------------------------
// LAGraph_Vector_TypeName: return the name of the GrB_Type of a GrB_Vector
//------------------------------------------------------------------------------

int LAGraph_Vector_TypeName
(
    // output:
    char *name,     // name of the type of the vector v (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    // input:
    GrB_Vector v,   // vector to query
    char *msg
)
{
    TYPENAME (v) ;
}

//------------------------------------------------------------------------------
// LAGraph_Scalar_TypeName: return the name of the GrB_Type of a GrB_Scalar
//------------------------------------------------------------------------------

int LAGraph_Scalar_TypeName
(
    // output:
    char *name,     // name of the type of the scalar s (user-provided array
                    // of size at least LAGRAPH_MAX_NAME_LEN).
    // input:
    GrB_Scalar s,   // scalar to query
    char *msg
)
{
    TYPENAME (s) ;
}

