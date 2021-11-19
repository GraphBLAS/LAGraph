//------------------------------------------------------------------------------
// LAGraph_Vector_isall: check two vectors
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// LAGraph_Vector_isall: check two vectors, contributed by Tim Davis, Texas A&M

// Applies a binary operator to two vectors A and B, and returns result = true
// if the pattern of A and B are identical, and if the result of C = A op B is
// true for all entries in C.

// See also LAGraph_isall for matrices.

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
    GrB_free (&C) ;

//****************************************************************************
GrB_Info LAGraph_Vector_isall      // return GrB_SUCCESS if successful
(
    bool *result,           // true if A == B, false if A != B or error
    GrB_Vector A,
    GrB_Vector B,
    GrB_BinaryOp op         // GrB_EQ_<type>, for the type of A and B,
                            // to check for equality.  Or use any desired
                            // operator.  The operator should return GrB_BOOL.
)
{
    GrB_Info info ;
    GrB_Vector C = NULL ;
    GrB_Index nrows1, nrows2, nvals, nvals1, nvals2 ;

    // check inputs
    if (result == NULL)
    {
        // error: required parameter, result, is NULL
        return (GrB_NULL_POINTER) ;
    }
    (*result) = false ;

    // check the size of A and B
    LAGRAPH_OK (GrB_Vector_size (&nrows1, A)) ;
    LAGRAPH_OK (GrB_Vector_size (&nrows2, B)) ;
    if (nrows1 != nrows2)
    {
        // # of rows differ
        // printf ("LAGraph_isall: rows differ\n") ;
        return (GrB_SUCCESS) ;
    }

    // check the # entries in A and B
    LAGRAPH_OK (GrB_Vector_nvals (&nvals1, A)) ;
    LAGRAPH_OK (GrB_Vector_nvals (&nvals2, B)) ;
    if (nvals1 != nvals2)
    {
        // # of entries differ
        // printf ("LAGraph_isall: nvals differ\n") ;
        return (GrB_SUCCESS) ;
    }

    // C = A .* B, where the pattern of C is the intersection of A and B
    LAGRAPH_OK (GrB_Vector_new (&C, GrB_BOOL, nrows1)) ;
    LAGRAPH_OK (GrB_eWiseMult (C, NULL, NULL, op, A, B, NULL)) ;

    // ensure C has the same number of entries as A and B
    LAGRAPH_OK (GrB_Vector_nvals (&nvals, C)) ;
    if (nvals != nvals1)
    {
        // pattern of A and B are different
        GrB_free (&C) ;
        return (GrB_SUCCESS) ;
    }

    // result = and (C)
    LAGRAPH_OK (GrB_reduce (result, NULL, GrB_LAND_MONOID_BOOL, C, NULL)) ;

    // printf ("isall : %d\n", *result) ;

    // free workspace and return result
    LAGraph_FREE_ALL;
    return (GrB_SUCCESS) ;
}
