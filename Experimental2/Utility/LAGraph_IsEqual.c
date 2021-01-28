//------------------------------------------------------------------------------
// LAGraph_IsEqual: check two matrices for exact equality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_isequal: check if two matrices are identically equal (same size,
// type, pattern, size, and values).

// For both methods, if the two matrices are GrB_FP32, GrB_FP64, GxB_FC32,
// or GxB_FC64 and have NaNs, then these functions will return false,
// since NaN == NaN is false.  To check for NaN equality (like
// isequalwithequalnans in MATLAB), use LAGraph_IsAll with a user-defined
// operator f(x,y) that returns true if x and y are both NaN.

#include "LG_internal.h"

int LAGraph_IsEqual         // returns 0 if successful, -1 if failure
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op,        // for A and B with arbitrary user-defined types.
                            // Ignored if A and B have built-in types.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Type atype, btype ;
    LG_CHECK (result == NULL, -1, "bad arguments") ;
    (*result) = false ;

    //--------------------------------------------------------------------------
    // compare the type of A and B
    //--------------------------------------------------------------------------

    GrB_TRY (GxB_Matrix_type (&atype, A)) ;
    GrB_TRY (GxB_Matrix_type (&btype, B)) ;
    if (atype != btype)
    {
        // types differ
        return (0) ;
    }

    //--------------------------------------------------------------------------
    // select the comparator operator
    //--------------------------------------------------------------------------


    GrB_BinaryOp compare ;
    // LAGraph_BinaryOp_Picker (&compare, "==", atype) ;
    // LAGraph_BinaryOp_Picker (&compare, "+", atype) ;
    // LAGraph_BinaryOp_Picker (&compare, "plus", atype) ;

    if      (atype == GrB_BOOL  ) compare = GrB_EQ_BOOL   ;
    else if (atype == GrB_INT8  ) compare = GrB_EQ_INT8   ;
    else if (atype == GrB_INT16 ) compare = GrB_EQ_INT16  ;
    else if (atype == GrB_INT32 ) compare = GrB_EQ_INT32  ;
    else if (atype == GrB_INT64 ) compare = GrB_EQ_INT64  ;
    else if (atype == GrB_UINT8 ) compare = GrB_EQ_UINT8  ;
    else if (atype == GrB_UINT16) compare = GrB_EQ_UINT16 ;
    else if (atype == GrB_UINT32) compare = GrB_EQ_UINT32 ;
    else if (atype == GrB_UINT64) compare = GrB_EQ_UINT64 ;
    else if (atype == GrB_FP32  ) compare = GrB_EQ_FP32   ;
    else if (atype == GrB_FP64  ) compare = GrB_EQ_FP64   ;
    else if (atype == GxB_FC32  ) compare = GxB_EQ_FC32   ;
    else if (atype == GxB_FC64  ) compare = GxB_EQ_FC64   ;
    else compare = op ;

    //--------------------------------------------------------------------------
    // compare the size, pattern, and values of A and B
    //--------------------------------------------------------------------------

    LAGraph_TRY (LAGraph_IsAll (result, A, B, compare, msg)) ;

    return (0) ;
}

