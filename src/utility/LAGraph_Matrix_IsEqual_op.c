//------------------------------------------------------------------------------
// LAGraph_Matrix_IsEqual_op: compare two matrices with a given op
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_Matrix_IsEqual_op: check if two matrices are equal (same
// size, structure, size, and values).

#define LG_FREE_WORK GrB_free (&C) ;

#include "LG_internal.h"

//------------------------------------------------------------------------------
// LAGraph_Matrix_IsEqual_op: compare two matrices using a given operator
//------------------------------------------------------------------------------

int LAGraph_Matrix_IsEqual_op
(
    // output:
    bool *result,           // true if A == B, false if A != B or error
    // input:
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op,        // comparator to use
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix C = NULL ;
    LG_ASSERT (op != NULL && result != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // check for NULL and aliased matrices
    //--------------------------------------------------------------------------

    if (A == NULL || B == NULL || A == B)
    {
        // two NULL matrices are identical, as are two aliased matrices
        (*result) = (A == B) ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // compare the size of A and B
    //--------------------------------------------------------------------------

    GrB_Index nrows1, ncols1, nrows2, ncols2 ;
    GrB_TRY (GrB_Matrix_nrows (&nrows1, A)) ;
    GrB_TRY (GrB_Matrix_nrows (&nrows2, B)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols1, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols2, B)) ;
    if (nrows1 != nrows2 || ncols1 != ncols2)
    {
        // dimensions differ
        (*result) = false ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // compare the # entries in A and B
    //--------------------------------------------------------------------------

    GrB_Index nvals1, nvals2 ;
    GrB_TRY (GrB_Matrix_nvals (&nvals1, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals2, B)) ;
    if (nvals1 != nvals2)
    {
        // # of entries differ
        (*result) = false ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // C = A .* B, where the structure of C is the intersection of A and B
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Matrix_new (&C, GrB_BOOL, nrows1, ncols1)) ;
    GrB_TRY (GrB_eWiseMult (C, NULL, NULL, op, A, B, NULL)) ;

    //--------------------------------------------------------------------------
    // ensure C has the same number of entries as A and B
    //--------------------------------------------------------------------------

    GrB_Index nvals ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, C)) ;
    if (nvals != nvals1)
    {
        // structure of A and B are different
        LG_FREE_WORK ;
        (*result) = false ;
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // result = and (C)
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, C, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

