//------------------------------------------------------------------------------
// LAGraph_pattern: return the pattern of a matrix (spones(A) in MATLAB)
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_pattern: return the pattern of a matrix (spones(A) in MATLAB)
// as a boolean matrix.

// SPEC: to do this in general for any user-defined types requires either (a)
// the user to create an operator z=f(x)=1, where z is boolean and x is the
// user type (LAGraph_TRUE_BOOL_Complex, for example), or (b)
// extractTuples(&I,&J,&X,A).  The latter requires X to be allocated of the
// right size, and then freed.  SuiteSparse allows X to be NULL but this is an
// extension to the spec. Determining the right size of X is difficult since
// there is no GrB_Type_size (see GxB_Type_size in SuiteSparse:GraphBLAS).

// As a result of these limitations, this method does not handle user-defined
// types, other than LAGraph_Complex (this function uses option (a) above).

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL \
    GrB_free (C) ;

GrB_Info LAGraph_pattern    // return GrB_SUCCESS if successful
(
    GrB_Matrix *C,          // a boolean matrix with the pattern of A
    GrB_Matrix A
)
{

    GrB_Info info ;
    GrB_Type type ;
    GrB_Index nrows, ncols ;

    // check inputs
    if (C == NULL)
    {
        // error: required parameter, result, is NULL
        return (GrB_NULL_POINTER) ;
    }
    (*C) = NULL ;

    // GxB_fprint (A, GxB_COMPLETE, stdout) ;

    // get the type and size of A
    LAGRAPH_OK (GxB_Matrix_type  (&type,  A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;

    // C = boolean matrix, the same size as A
    LAGRAPH_OK (GrB_Matrix_new (C, GrB_BOOL, nrows, ncols)) ;

    // C = spones (A), typecasting to bool
    if (type == LAGraph_Complex)
    {
        // the LAGraph_TRUE_BOOL_Complex operator returns boolean true,
        // and has an input of type LAGraph_Complex (which it ignores).
        LAGRAPH_OK (GrB_apply (*C, NULL, NULL, LAGraph_TRUE_BOOL_Complex,
            A, NULL)) ;
    }
    else
    {
        // this works for all built-in types, which are first typecasted to
        // boolean ... and then ignored by the operator anyway ...
        LAGRAPH_OK (GrB_apply (*C, NULL, NULL, LAGraph_TRUE_BOOL, A, NULL)) ;
    }

    // free workspace and return result
    return (GrB_SUCCESS) ;
}

