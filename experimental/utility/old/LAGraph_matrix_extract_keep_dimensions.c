//------------------------------------------------------------------------------
// LAGraph_matrix_extract_keep_dimensions: extract submatrix but keep the
// dimensions of the original matrix
// ------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// LAGraph_Matrix_extract_keep_dimensions: Contributed by Gabor Szarnyas.
// Budapest University of Technology and Economics
// (with accented characters: G\'{a}bor Sz\'{a}rnyas).

// Compute the

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL            \
{                                   \
    GrB_free(&C) ;                  \
    GrB_free(&type) ;               \
}

//****************************************************************************
typedef struct
{
    const GrB_Index nv; // number of vertices
    const bool* Vdense; // array denoting whether a vertex should be kept
} Vdense_struct_type;


bool select_submatrix_elements_fun(
const GrB_Index i, const GrB_Index j,
const void *x, const void *thunk)
{
    Vdense_struct_type* indices = (Vdense_struct_type*) (thunk);
    return indices->Vdense[i] && indices->Vdense[j];
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_Matrix_extract_keep_dimensions // extract submatrix but keep
                                                // the dimensions of the
                                                // original matrix
(
    GrB_Matrix *Chandle,         // output matrix
    const GrB_Matrix A,          // input matrix
    const GrB_Index *Vsparse,    // sorted list of vertex indices
    const bool *Vdense,          // boolean array of vertices
    GrB_Index nv                 // number of vertex indices
)
{
#if !(LG_SUITESPARSE)
    return GrB_PANIC;
#else
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Type type;
    GrB_Index n ;
    GrB_Matrix C = NULL ;

    LAGRAPH_OK (GxB_Matrix_type(&type, A));
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A));
    LAGRAPH_OK (GrB_Matrix_new (&C, type, n, n));

    if (Vsparse == NULL && Vdense == NULL)
    {
        LAGRAPH_ERROR("Both Vsparse and Vdense are set to NULL", GrB_NULL_POINTER)
    }

    if (Vsparse == NULL) // use Vdense and GxB_select
    {
        Vdense_struct_type vdense_struct = {.nv = nv, .Vdense = Vdense};

        GrB_Type Vdense_type;
        LAGRAPH_OK (GrB_Type_new(&Vdense_type, sizeof(vdense_struct)));

        GxB_Scalar vdense_thunk;
        LAGRAPH_OK (GxB_Scalar_new(&vdense_thunk, Vdense_type));
        LAGRAPH_OK (GxB_Scalar_setElement(vdense_thunk, (void*) &vdense_struct));

        GxB_SelectOp select_submatrix_elements_op;
        LAGRAPH_OK (GxB_SelectOp_new(
                        &select_submatrix_elements_op,
                        select_submatrix_elements_fun,
                        NULL, Vdense_type));
        LAGRAPH_OK (GxB_select(C, NULL, NULL, select_submatrix_elements_op,
                               A, vdense_thunk, NULL));

        GrB_free(&select_submatrix_elements_op); select_submatrix_elements_op = NULL;
        GrB_free(&vdense_thunk); vdense_thunk = NULL;
        GrB_free(&Vdense_type); Vdense_type = NULL;
    }
    else
    {
        GrB_Matrix D; // diagonal matrix used to select rows/columns

        LAGRAPH_OK (GrB_Matrix_new(&D, GrB_BOOL, n, n));

        bool* X = LAGraph_Malloc(nv, sizeof(GrB_BOOL)) ;
        if (X == NULL)
        {
            LAGRAPH_ERROR("out of memory", GrB_OUT_OF_MEMORY)
        }
        int nthreads;
        LAGraph_GetNumThreads(&nthreads, NULL) ;
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (GrB_Index i = 0; i < nv; i++)
        {
            X[i] = true;
        }
        LAGRAPH_OK (GrB_Matrix_build(D, Vsparse, Vsparse, X, nv, GrB_LOR));

        GxB_Format_Value A_format;
        LAGRAPH_OK(GxB_get(A, GxB_FORMAT, &A_format))
        if (A_format == GxB_BY_ROW) // C = (D*A)*D
        {
            LAGRAPH_OK (GrB_mxm(C, NULL, NULL, GxB_ANY_SECOND_FP64, D, A, NULL));
            LAGRAPH_OK (GrB_mxm(C, NULL, NULL, GxB_ANY_FIRST_FP64, C, D, NULL));
        }
        else // A_format == GxB_BY_COL: C = D*(A*D)
        {
            LAGRAPH_OK (GrB_mxm(C, NULL, NULL, GxB_ANY_FIRST_FP64, A, D, NULL));
            LAGRAPH_OK (GrB_mxm(C, NULL, NULL, GxB_ANY_SECOND_FP64, D, C, NULL));
        }

        GrB_free(&D); D = NULL;
    }
    (*Chandle) = C ;

    return (GrB_SUCCESS) ;
#endif
}
