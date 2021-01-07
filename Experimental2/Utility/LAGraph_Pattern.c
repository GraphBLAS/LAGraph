//------------------------------------------------------------------------------
// LAGraph_Pattern: return the pattern of a matrix (spones(A) in MATLAB)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// LAGraph_pattern: return the pattern of a matrix (spones(A) in MATLAB)
// as a boolean matrix. Contributed by Tim Davis and Scott Kolodziej, Texas A&M.

#include "LAGraph_Internal.h"

#define LAGRAPH_FREE_ALL \
    GrB_free (C) ;

int LAGraph_Pattern     // return 0 if successful, -1 if failure
(
    GrB_Matrix *C,      // a boolean matrix with the pattern of A
    GrB_Matrix A,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LAGraph_CLEAR_MSG ;
    GrB_Index nrows, ncols ;
    LAGraph_CHECK (C == NULL, -1, "&C is NULL") ;
    LAGraph_CHECK (A == NULL, -1, "A is NULL") ;
    (*C) = NULL ;

    // GxB_fprint (A, GxB_COMPLETE, stdout) ;

    //--------------------------------------------------------------------------
    // get the size of A
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

    //--------------------------------------------------------------------------
    // C = spones (A), typecasting to bool
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_new (C, GrB_BOOL, nrows, ncols)) ;
    GrB_TRY (GrB_apply (*C, NULL, NULL, GxB_ONE_BOOL, A, NULL)) ;

    // free workspace and return result
    return (0) ;
}

