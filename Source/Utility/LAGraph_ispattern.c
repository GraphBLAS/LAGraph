//------------------------------------------------------------------------------
// LAGraph_ispattern: check if a matrix is all 1
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_ispattern: check if a matrix values are all equal to 1.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
    GrB_free (&C) ;

GrB_Info LAGraph_ispattern  // return GrB_SUCCESS if successful
(
    bool *result,           // true if A is all one, false otherwise
    GrB_Matrix A,
    GrB_UnaryOp userop      // for A with arbitrary user-defined type.
                            // Ignored if A and B are of built-in types or
                            // LAGraph_Complex.
)
{

    GrB_Matrix C = NULL ;
    GrB_Type type ;
    GrB_Index nrows, ncols ;

    // check inputs
    if (result == NULL)
    {
        // error: required parameter, result, is NULL
        return (GrB_NULL_POINTER) ;
    }
    (*result) = false ;

    // GxB_fprint (A, GxB_COMPLETE, stdout) ;

    // get the type and size of A
    LAGRAPH_OK (GxB_Matrix_type  (&type,  A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;

    if (type == GrB_BOOL)
    {
        // result = and (A)
        LAGRAPH_OK (GrB_reduce (result, NULL, LAGraph_LAND_MONOID, A, NULL)) ;
    }
    else
    {

        // select the unary operator
        GrB_UnaryOp op = NULL ;
        if      (type == GrB_INT8  ) op = LAGraph_ISONE_INT8 ;
        else if (type == GrB_INT16 ) op = LAGraph_ISONE_INT16 ;
        else if (type == GrB_INT32 ) op = LAGraph_ISONE_INT32 ;
        else if (type == GrB_INT64 ) op = LAGraph_ISONE_INT64 ;
        else if (type == GrB_UINT8 ) op = LAGraph_ISONE_UINT8 ;
        else if (type == GrB_UINT16) op = LAGraph_ISONE_UINT16 ;
        else if (type == GrB_UINT32) op = LAGraph_ISONE_UINT32 ;
        else if (type == GrB_UINT64) op = LAGraph_ISONE_UINT64 ;
        else if (type == GrB_FP32  ) op = LAGraph_ISONE_FP32   ;
        else if (type == GrB_FP64  ) op = LAGraph_ISONE_FP64   ;
        else if (type == LAGraph_Complex) op = LAGraph_ISONE_Complex   ;
        else op = userop ;

        if (op == NULL)
        {
            printf ("LAGraph_ispattern: userop is NULL\n") ;
            return (GrB_NULL_POINTER) ;
        }

        // C = isone (A)
        LAGRAPH_OK (GrB_Matrix_new (&C, GrB_BOOL, nrows, ncols)) ;

        LAGRAPH_OK (GrB_apply (C, NULL, NULL, op, A, NULL)) ;
        // GxB_fprint (C, GxB_COMPLETE, stdout) ;

        // result = and (C)
        LAGRAPH_OK (GrB_reduce (result, NULL, LAGraph_LAND_MONOID, C, NULL)) ;
    }

    // free workspace and return result
    GrB_free (&C) ;
    return (GrB_SUCCESS) ;
}

