//------------------------------------------------------------------------------
// LAGraph_IsEqual: check two matrices for exact equality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_IsEqual: check if two matrices are identically equal (same size,
// type, structure, size, and values).

// If the two matrices are GrB_FP32, GrB_FP64, GxB_FC32, or GxB_FC64 and have
// NaNs, then these functions will return false, since NaN == NaN is false.

#define LAGraph_FREE_WORK GrB_free (&C) ;

#include "LG_internal.h"

//------------------------------------------------------------------------------
// LAGraph_IsEqual_type: compare using GrB_EQ_type operator
//------------------------------------------------------------------------------

int LAGraph_IsEqual_type    // returns 0 if successful, < 0 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_Type type,          // use GrB_EQ_type operator to compare A and B
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix C = NULL ;
    LG_CHECK (type == NULL || result == NULL, -1001, "required input is NULL") ;

    //--------------------------------------------------------------------------
    // check for NULL and aliased matrices
    //--------------------------------------------------------------------------

    if (A == NULL || B == NULL || A == B)
    {
        // two NULL matrices are identical, as are two aliased matrices
        (*result) = (A == B) ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // compare the size of A and B
    //--------------------------------------------------------------------------

    GrB_Index nrows1, ncols1, nrows2, ncols2 ;
    GrB_TRY (GrB_Matrix_nrows (&nrows1, A)) ;
    GrB_TRY (GrB_Matrix_nrows (&nrows2, B)) ;
    if (nrows1 != nrows2)
    {
        // # of rows differ
        (*result) = false ;
        return (0) ;
    }

    GrB_TRY (GrB_Matrix_ncols (&ncols1, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols2, B)) ;
    if (ncols1 != ncols2)
    {
        // # of cols differ
        (*result) = false ;
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
        (*result) = false ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // get the GrB_EQ_type operator
    //--------------------------------------------------------------------------

    GrB_BinaryOp op ;
    // LAGraph_BinaryOp_Picker (&op, "==", type) ;
    // LAGraph_BinaryOp_Picker (&op, "+", type) ;
    // LAGraph_BinaryOp_Picker (&op, "plus", type) ;
    if      (type == GrB_BOOL  ) op = GrB_EQ_BOOL   ;
    else if (type == GrB_INT8  ) op = GrB_EQ_INT8   ;
    else if (type == GrB_INT16 ) op = GrB_EQ_INT16  ;
    else if (type == GrB_INT32 ) op = GrB_EQ_INT32  ;
    else if (type == GrB_INT64 ) op = GrB_EQ_INT64  ;
    else if (type == GrB_UINT8 ) op = GrB_EQ_UINT8  ;
    else if (type == GrB_UINT16) op = GrB_EQ_UINT16 ;
    else if (type == GrB_UINT32) op = GrB_EQ_UINT32 ;
    else if (type == GrB_UINT64) op = GrB_EQ_UINT64 ;
    else if (type == GrB_FP32  ) op = GrB_EQ_FP32   ;
    else if (type == GrB_FP64  ) op = GrB_EQ_FP64   ;
    #if 0
    else if (type == GxB_FC32  ) op = GxB_EQ_FC32   ;
    else if (type == GxB_FC64  ) op = GxB_EQ_FC64   ;
    #endif
    else
    {
        LG_CHECK (true, -1002, "unsupported type") ;
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
        LAGraph_FREE_WORK ;
        (*result) = false ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // result = and (C)
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, C, NULL)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_IsEqual: compare using GrB_EQ_type operator; auto type selection
//------------------------------------------------------------------------------

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK ;

int LAGraph_IsEqual         // returns 0 if successful, < 0 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (A == NULL, -1001, "required input is NULL") ;

    //--------------------------------------------------------------------------
    // determine the type
    //--------------------------------------------------------------------------

    GrB_Type type ;
    #if LG_SUITESPARSE
        // SuiteSparse:GraphBLAS: query the type and compare accordingly
        GrB_TRY (GxB_Matrix_type (&type, A)) ;
    #else
        // no way to determine the type with pure GrB*; compare as if FP64
        type = GrB_FP64 ;
    #endif

    //--------------------------------------------------------------------------
    // compare A and B
    //--------------------------------------------------------------------------

    return (LAGraph_IsEqual_type (result, A, B, type, msg)) ;
}

