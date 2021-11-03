//------------------------------------------------------------------------------
// LAGraph_Vector_wait: interface to GrB_Vector_wait
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Vector_wait     // wait on a vector
(
    GrB_Vector v,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (v == NULL, -1001, "A is NULL") ;

    //--------------------------------------------------------------------------
    // wait on the vector
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_wait (v, GrB_MATERIALIZE)) ;
    return (0) ;
}
