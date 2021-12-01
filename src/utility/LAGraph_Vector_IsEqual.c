//------------------------------------------------------------------------------
// LAGraph_Vector_IsEqual: check two vectors for exact equality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_Vector_IsEqual, contributed by Tim Davis, Texas A&M

// Checks if two vectors are identically equal (same size,
// type, pattern, size, and values).  Checking for the same type requires the
// GxB_Vector_type function, which is an extension in SuiteSparse:GraphBLAS.
// For the standard API, there is no way to determine the type of a vector.

// See also LAGraph_IsEqual.

// For both methods, if the two vectors are GrB_FP32, GrB_FP64, or related,
// and have NaNs, then these functions will return false,
// since NaN == NaN is false.  To check for NaN equality (like
// isequalwithequalnans in MATLAB), use LAGraph_isall with a user-defined
// operator f(x,y) that returns true if x and y are both NaN.

#define LAGraph_FREE_WORK GrB_free (&C) ;

#include "LG_internal.h"
#include <LAGraphX.h>

//****************************************************************************
//****************************************************************************

GrB_Info LAGraph_Vector_IsEqual_op    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp userop,    // comparator to use (required)
    char *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector C = NULL ;
    LG_CHECK (userop == NULL || result == NULL, -1001,
        "required input is NULL") ;

    GrB_Info info ;

    //--------------------------------------------------------------------------
    // check for NULL and aliased vectors
    //--------------------------------------------------------------------------

    if (A == NULL || B == NULL || A == B)
    {
        // two NULL vectors are identical, as are two aliased matrices
        (*result) = (A == B) ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // compare the size of A and B
    //--------------------------------------------------------------------------

    GrB_Index nrows1, nrows2;
    GrB_TRY (GrB_Vector_size (&nrows1, A)) ;
    GrB_TRY (GrB_Vector_size (&nrows2, B)) ;
    if (nrows1 != nrows2)
    {
        // # of rows differ
        (*result) = false ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // compare the # entries in A and B
    //--------------------------------------------------------------------------

    GrB_Index nvals1, nvals2 ;
    GrB_TRY (GrB_Vector_nvals (&nvals1, A)) ;
    GrB_TRY (GrB_Vector_nvals (&nvals2, B)) ;
    if (nvals1 != nvals2)
    {
        // # of entries differ
        (*result) = false ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // C = A .* B, where the pattern of C is the intersection of A and B
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&C, GrB_BOOL, nrows1)) ;
    GrB_TRY (GrB_eWiseMult (C, NULL, NULL, userop, A, B, NULL)) ;

    //--------------------------------------------------------------------------
    // ensure C has the same number of entries as A and B
    //--------------------------------------------------------------------------

    GrB_Index nvals ;
    GrB_TRY (GrB_Vector_nvals (&nvals, C)) ;
    if (nvals != nvals1)
    {
        // pattern of A and B are different
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

    // return result
    LAGraph_FREE_WORK;
    return 0;
}


//****************************************************************************
//****************************************************************************
GrB_Info LAGraph_Vector_IsEqual_type    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_Type   type,         // use GrB_EQ_type operator to compare A and B
    char *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector C = NULL ;
    LG_CHECK (type == NULL || result == NULL, -1001, "required input is NULL") ;

    GrB_Info info ;

    //--------------------------------------------------------------------------
    // check for NULL and aliased vectors
    //--------------------------------------------------------------------------

    if (A == NULL || B == NULL || A == B)
    {
        // two NULL vectors are identical, as are two aliased matrices
        (*result) = (A == B) ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // compare the size of A and B
    //--------------------------------------------------------------------------

    GrB_Index nrows1, nrows2;
    GrB_TRY (GrB_Vector_size (&nrows1, A)) ;
    GrB_TRY (GrB_Vector_size (&nrows2, B)) ;
    if (nrows1 != nrows2)
    {
        // # of rows differ
        (*result) = false ;
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // compare the # entries in A and B
    //--------------------------------------------------------------------------

    GrB_Index nvals1, nvals2 ;
    GrB_TRY (GrB_Vector_nvals (&nvals1, A)) ;
    GrB_TRY (GrB_Vector_nvals (&nvals2, B)) ;
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
    // select the comparator operator
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
    // C = A .* B, where the pattern of C is the intersection of A and B
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&C, GrB_BOOL, nrows1)) ;
    GrB_TRY (GrB_eWiseMult (C, NULL, NULL, op, A, B, NULL)) ;

    //--------------------------------------------------------------------------
    // ensure C has the same number of entries as A and B
    //--------------------------------------------------------------------------

    GrB_Index nvals ;
    GrB_TRY (GrB_Vector_nvals (&nvals, C)) ;
    if (nvals != nvals1)
    {
        // pattern of A and B are different
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

    // return result
    LAGraph_FREE_WORK;
    return 0;
}

//------------------------------------------------------------------------------
// LAGraph_IsEqual: compare using GrB_EQ_type operator; auto type selection
//------------------------------------------------------------------------------

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK ;

int LAGraph_Vector_IsEqual         // returns 0 if successful, < 0 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
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
        GrB_TRY (GxB_Vector_type (&type, A)) ;
    #else
        // no way to determine the type with pure GrB*; compare as if FP64
        type = GrB_FP64 ;
    #endif

    //--------------------------------------------------------------------------
    // compare A and B
    //--------------------------------------------------------------------------

    return (LAGraph_Vector_IsEqual_type (result, A, B, type, msg)) ;
}
