//------------------------------------------------------------------------------
// LAGraph_SizeOfType: return the sizeof(...) of a GraphBLAS GrB_Type
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

#include "LG_internal.h"

int LAGraph_SizeOfType
(
    // output:
    size_t *size,   // size of the type
    // input:
    GrB_Type type,  // GraphBLAS type
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (type != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (size != NULL, GrB_NULL_POINTER) ;
    (*size) = 0 ;

    //--------------------------------------------------------------------------
    // determine the size of the type
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Type_get_SIZE (type, size, GrB_SIZE)) ;
    return (GrB_SUCCESS) ;
}

