//------------------------------------------------------------------------------
// LAGraph_Matrix_wait: interface to GrB_Matrix_wait
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Matrix_wait     // wait on a matrix
(
    GrB_Matrix A,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (A == NULL, -1001, "A is NULL") ;

    //--------------------------------------------------------------------------
    // wait on the matrix
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_wait (A, GrB_MATERIALIZE)) ;
    return (0) ;
}

