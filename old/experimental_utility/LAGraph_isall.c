//------------------------------------------------------------------------------
// LAGraph_isall: check two matrices
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// LAGraph_isall: check two matrices.  Contributed by Tim Davis, Texas A&M

// Applies a binary operator to two matrices A and B, and returns result = true
// if the pattern of A and B are identical, and if the result of C = A op B is
// true for all entries in C.

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
    GrB_free (&C) ;

//****************************************************************************
GrB_Info LAGraph_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Matrix A,
    GrB_Matrix B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
)
{
    GrB_Info info ;
    GrB_Matrix C = NULL ;
    GrB_Index nrows1, ncols1, nrows2, ncols2, nvals, nvals1, nvals2 ;

    // check inputs
    if (result == NULL)
    {
        // error: required parameter, result, is NULL
        // printf ("LAGraph_isall: bad args \n") ;
        return (GrB_NULL_POINTER) ;
    }
    (*result) = false ;

    // check the size of A and B
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows1, A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows2, B)) ;
    if (nrows1 != nrows2)
    {
        // # of rows differ
        return (GrB_SUCCESS) ;
    }

    LAGRAPH_OK (GrB_Matrix_ncols (&ncols1, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols2, B)) ;
    if (ncols1 != ncols2)
    {
        // # of cols differ
        return (GrB_SUCCESS) ;
    }

    // check the # entries in A and B
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals1, A)) ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals2, B)) ;
    if (nvals1 != nvals2)
    {
        // # of entries differ
        return (GrB_SUCCESS) ;
    }

    // C = A .* B, where the pattern of C is the intersection of A and B
    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_BOOL, nrows1, ncols1)) ;
    LAGRAPH_OK (GrB_eWiseMult (C, NULL, NULL, op, A, B, NULL)) ;

    // ensure C has the same number of entries as A and B
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, C)) ;
    if (nvals != nvals1)
    {
        // pattern of A and B are different
        GrB_free (&C) ;
        return (GrB_SUCCESS) ;
    }

    // result = and (C)
    LAGRAPH_OK (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, C, NULL)) ;

    // free workspace and return result
    GrB_free (&C) ;
    return (GrB_SUCCESS) ;
}
