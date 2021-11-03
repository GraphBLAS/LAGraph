//------------------------------------------------------------------------------
// LAGraph_Structure: return the structure of a matrix (spones(A) in MATLAB)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis and Scott Kolodziej, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_Structure: return the structure of a matrix as a boolean matrix,
// where C(i,j)=true if the entry A(i,j) is present in the matrix A.

#define LAGraph_FREE_ALL GrB_free (C) ;

#include "LG_internal.h"

int LAGraph_Structure   // return 0 if successful, -1 if failure
(
    GrB_Matrix *C,      // a boolean matrix with same structure of A
    GrB_Matrix A,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Index nrows, ncols ;
    LG_CHECK (C == NULL, -1, "&C is NULL") ;
    LG_CHECK (A == NULL, -1, "A is NULL") ;
    (*C) = NULL ;

    //--------------------------------------------------------------------------
    // get the size of A
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;

    //--------------------------------------------------------------------------
    // C<s(A)> = true
    // C{A} = true
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_new (C, GrB_BOOL, nrows, ncols)) ;
    GrB_TRY (GrB_assign (*C, A, NULL, (bool) true,
        GrB_ALL, nrows, GrB_ALL, nrows, GrB_DESC_S)) ;

    return (0) ;
}

