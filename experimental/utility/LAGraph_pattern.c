//------------------------------------------------------------------------------
// LAGraph_pattern: return the pattern of a matrix (spones(A) in MATLAB)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite
// FIXME: remove this anduse LAGraph_Structure.

// LAGraph_pattern: return the pattern of a matrix (spones(A) in MATLAB)
// as a boolean matrix. Contributed by Tim Davis and Scott Kolodziej, Texas A&M.

// SPEC: to do this in general for any user-defined types requires either (a)
// the user to create an operator z=f(x)=1, where z is boolean and x is the
// user type, or (b) extractTuples(&I,&J,&X,A).  The latter requires X to be
// allocated of the right size, and then freed.  SuiteSparse allows X to be
// NULL but this is an extension to the spec. Determining the right size of X
// is difficult since there is no GrB_Type_size (see GxB_Type_size in
// SuiteSparse:GraphBLAS).

// As a result of these limitations, this method does not handle user-defined
// types.

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL               \
{                                      \
    GrB_free (C);                      \
    GrB_free (&LAGraph_TRUE_BOOL);     \
}

//****************************************************************************
// unary operators that return boolean true
#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

void LAGraph_true_bool
(
    bool *z,
    const bool *x       // ignored
)
{
    (*z) = true ;
}

//****************************************************************************
GrB_Info LAGraph_pattern    // return GrB_SUCCESS if successful
(
    GrB_Matrix *C,          // a boolean matrix with the pattern of A
    GrB_Matrix A,
    GrB_Type C_type
)
{
    GrB_Info info;
    GrB_Type type;
    GrB_Index nrows, ncols;
    GrB_UnaryOp LAGraph_TRUE_BOOL = NULL ;

    // check inputs
    if (C == NULL)
    {
        // error: required parameter, result, is NULL
        return (GrB_NULL_POINTER);
    }
    (*C) = NULL ;

    if (C_type == GrB_NULL)
    {
        C_type = GrB_BOOL;
    }

    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_TRUE_BOOL,
                                 F_UNARY (LAGraph_true_bool),
                                 GrB_BOOL, GrB_BOOL)) ;

    // GxB_fprint (A, GxB_COMPLETE, stdout) ;

    // get the type and size of A
    //LAGRAPH_OK (GxB_Matrix_type  (&type,  A));
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A));

    // C = boolean matrix, the same size as A
    // T = GrB_BOOL by default
    LAGRAPH_OK (GrB_Matrix_new (C, C_type, nrows, ncols));

    // this works for all built-in types, which are first typecasted to
    // boolean ... and then ignored by the operator anyway ...
    LAGRAPH_OK (GrB_apply(*C, NULL, NULL, LAGraph_TRUE_BOOL, A, NULL));

    // free workspace and return result
    return (GrB_SUCCESS);
}
