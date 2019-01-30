//------------------------------------------------------------------------------
// LAGraph_isequal: check two matrices for exact equality
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_isequal: check if two matrices are identically equal (same size,
// type, pattern, size, and values).  Checking for the same type requires the
// GxB_Matrix_type function, which is an extension in SuiteSparse:GraphBLAS.
// For the standard API, there is no way to determine the type of a matrix.

// TODO: add GrB_Matrix_Type to the GraphBLAS spec.

// For both methods, if the two matrices are GrB_FP32, GrB_FP64, or
// LAGraph_Complex, and have NaNs, then these functions will return false,
// since NaN == NaN is false.  To check for NaN equality (like
// isequalwithequalnans in MATLAB), use LAGraph_isall with a user-defined
// operator f(x,y) that returns true if x and y are both NaN.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL ;

GrB_Info LAGraph_isequal    // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp userop     // for A and B with arbitrary user-defined types.
                            // Ignored if A and B are of built-in types or
                            // LAGraph_Complex.
)
{

    GrB_Type atype, btype ;
    GrB_BinaryOp op ;

    // check inputs
    if (result == NULL)
    {
        // error: required parameter, result, is NULL
        printf ("LAGraph_isequal: bad args\n") ;
        return (GrB_NULL_POINTER) ;
    }
    (*result) = false ;

    // check the type of A and B
    LAGRAPH_OK (GxB_Matrix_type (&atype, A)) ;
    LAGRAPH_OK (GxB_Matrix_type (&btype, B)) ;
    if (atype != btype)
    {
        // types differ
        // printf ("LAGraph_isequal: types differ\n") ;
        return (GrB_SUCCESS) ;
    }

    // select the comparator operator
    if      (atype == GrB_BOOL  ) op = GrB_EQ_BOOL ;
    else if (atype == GrB_INT8  ) op = GrB_EQ_INT8 ;
    else if (atype == GrB_INT16 ) op = GrB_EQ_INT16 ;
    else if (atype == GrB_INT32 ) op = GrB_EQ_INT32 ;
    else if (atype == GrB_INT64 ) op = GrB_EQ_INT64 ;
    else if (atype == GrB_UINT8 ) op = GrB_EQ_UINT8 ;
    else if (atype == GrB_UINT16) op = GrB_EQ_UINT16 ;
    else if (atype == GrB_UINT32) op = GrB_EQ_UINT32 ;
    else if (atype == GrB_UINT64) op = GrB_EQ_UINT64 ;
    else if (atype == GrB_FP32  ) op = GrB_EQ_FP32   ;
    else if (atype == GrB_FP64  ) op = GrB_EQ_FP64   ;
    else if (atype == LAGraph_Complex) op = LAGraph_EQ_Complex   ;
    else op = userop ;

    // check the size, pattern, and values of A and B
    LAGRAPH_OK (LAGraph_isall (result, A, B, op)) ;

    // printf ("result: %d\n", *result) ;

    // return result
    return (GrB_SUCCESS) ;
}

