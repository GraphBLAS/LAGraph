//------------------------------------------------------------------------------
// LAGraph_IsAll: compare two matrices
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// Applies a binary operator to two matrices A and B, and returns result = true
// if the pattern of A and B are identical, and if the result of C = A op B is
// true for all entries in C.

#include "LG_internal.h"

#define LAGRAPH_FREE_WORK   \
    GrB_free (&C) ;         \

int LAGraph_IsAll           // returns 0 if successful, -1 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op,        // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Matrix C = NULL ;
    LG_CHECK (result == NULL, -1, "bad args") ;
    (*result) = false ;

    //--------------------------------------------------------------------------
    // compare the size of A and B
    //--------------------------------------------------------------------------

    GrB_Index nrows1, ncols1, nrows2, ncols2 ;
    GrB_TRY (GrB_Matrix_nrows (&nrows1, A)) ;
    GrB_TRY (GrB_Matrix_nrows (&nrows2, B)) ;
    if (nrows1 != nrows2)
    {
        // # of rows differ
        return (0) ;
    }

    GrB_TRY (GrB_Matrix_ncols (&ncols1, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols2, B)) ;
    if (ncols1 != ncols2)
    {
        // # of cols differ
        return (0) ;
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
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // C = A .* B, where the pattern of C is the intersection of A and B
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
        // pattern of A and B are different
        LAGRAPH_FREE_WORK ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // result = and (C)
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, C, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_WORK ;
    return (0) ;
}

