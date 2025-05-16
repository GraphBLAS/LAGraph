//----------------------------------------------------------------------------
// LAGraph/experimental/test/LG_check_argminmax.c: simple arg min/max method
// ----------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Tim Davis, Texas A&M University

// ----------------------------------------------------------------------------

#include "LG_Xtest.h"
#include "LG_internal.h"

#undef  LG_FREE_WORK
#define LG_FREE_WORK                    \
{                                       \
    LAGraph_Free ((void **) &I, msg) ;  \
    LAGraph_Free ((void **) &J, msg) ;  \
    LAGraph_Free ((void **) &X, msg) ;  \
    LAGraph_Free ((void **) &v, msg) ;  \
    LAGraph_Free ((void **) &p, msg) ;  \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                     \
{                                       \
    LG_FREE_WORK ;                      \
    GrB_free (x_result) ;               \
    GrB_free (p_result) ;               \
}

//------------------------------------------------------------------------------
// LG_check_argminmax: compute argmin/max of each row/column of A
//------------------------------------------------------------------------------

// This method is single threaded and thus slow (like many LG_check_* methods).
// For simplicity, this method does all its computations in FP64, and then
// typecasts its results when done.  This can cause loss of precision if the
// input matrix is of type int64 or uint64 and has entries larger than about
// 2^52.  It can compute a different result for the position p if dim==0;
// otherwise, its results should match LAGraph_argminmax.

int LG_check_argminmax
(
    // output
    GrB_Vector *x_result,       // min/max value in each row/col of A
    GrB_Vector *p_result,       // index of min/max value in each row/col of A
    // input
    GrB_Matrix A,
    int dim,                    // dim=1: cols of A, dim=2: rows of A
                                // dim=0: return a scalar A(i,j) in x(0) and
                                // its row and col in p.
    bool is_min,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    uint64_t *I = NULL, *J = NULL, *p = NULL ;
    double *X = NULL, *v = NULL ;

    LG_ASSERT (A != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (x_result != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (p_result != NULL, GrB_NULL_POINTER) ;

    (*x_result) = NULL ;
    (*p_result) = NULL ;

    //--------------------------------------------------------------------------
    // extract the entries from the matrix
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols, nvals, n, np ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;

    switch (dim)
    {
        case 0: n = 1     ; np = 2 ; break ;
        case 1: n = ncols ; np = n ; break ;
        case 2: n = nrows ; np = n ; break ;
        default: LG_ASSERT (false, GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // extract the entries from the matrix
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Malloc ((void **) &I, nvals, sizeof (uint64_t), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &J, nvals, sizeof (uint64_t), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &X, nvals, sizeof (double), msg)) ;

    GRB_TRY (GrB_Matrix_extractTuples_FP64 (I, J, X, &nvals, A)) ;

    GrB_Type type = NULL ;
    #if LAGRAPH_SUITSPARSE
    GRB_TRY (GxB_Matrix_type (&type, A)) ;
    #else
    char typename [LAGRAPH_MAX_NAME_LEN+1] ;
    LG_TRY (LAGraph_Matrix_TypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;
    #endif

    //--------------------------------------------------------------------------
    // allocate temporary arrays for the results
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Calloc ((void **) &v, n, sizeof (double), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &p, np, sizeof (uint64_t), msg)) ;
    for (int64_t k = 0 ; k < np ; k++)
    {
        // kth entry is not yet present
        p [k] = -1 ;
    }

    //--------------------------------------------------------------------------
    // compute the argmin/argmax
    //--------------------------------------------------------------------------

    if (dim == 0)
    {

        //----------------------------------------------------------------------
        // global argmin/argmax
        //----------------------------------------------------------------------

        if (nvals > 0)
        {
            // (v,p) = first entry in the matrix
            v [0] = X [0] ;
            p [0] = I [0] ;
            p [1] = J [0] ;

            for (int64_t k = 1 ; k < nvals ; k++)
            {
                if (is_min)
                {
                    // argmin
                    if (X [k] < v [0])
                    {
                        v [0] = X [k] ;
                        p [0] = I [k] ;
                        p [1] = J [k] ;
                    }
                }
                else
                {
                    // argmax
                    if (X [k] > v [0])
                    {
                        v [0] = X [k] ;
                        p [0] = I [k] ;
                        p [1] = J [k] ;
                    }
                }
            }

        }

    }
    else if (dim == 1)
    {

        //----------------------------------------------------------------------
        // argmin/argmax of each column of A
        //----------------------------------------------------------------------

        for (int64_t k = 0 ; k < nvals ; k++)
        {
            int64_t i = I [k] ;
            int64_t j = J [k] ;
            double x = X [k] ;
            if (p [j] == -1)
            {
                // first entry seen in column j
                v [j] = x ;
                p [j] = i ;
            }
            else if (is_min)
            {
                if (x < v [j])
                {
                    v [j] = x ;
                    p [j] = i ;
                }
            }
            else
            {
                if (x > v [j])
                {
                    v [j] = x ;
                    p [j] = i ;
                }
            }
        }

    }
    else // (dim == 2)
    {

        //----------------------------------------------------------------------
        // argmin/argmax of each row of A
        //----------------------------------------------------------------------

        for (int64_t k = 0 ; k < nvals ; k++)
        {
            int64_t i = I [k] ;
            int64_t j = J [k] ;
            double x = X [k] ;
            if (p [i] == -1)
            {
                // first entry seen in row i
                v [i] = x ;
                p [i] = j ;
            }
            else if (is_min)
            {
                if (x < v [i])
                {
                    v [i] = x ;
                    p [i] = j ;
                }
            }
            else
            {
                if (x > v [i])
                {
                    v [i] = x ;
                    p [i] = j ;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // create the output vectors and fill with the results
    //--------------------------------------------------------------------------

    LG_TRY (GrB_Vector_new (x_result, type, n)) ;
    LG_TRY (GrB_Vector_new (p_result, GrB_INT64, np)) ;

    if (dim == 0)
    {
        if (nvals > 0)
        {
            GRB_TRY (GrB_Vector_setElement_FP64  (*x_result, v [0], 0)) ;
            GRB_TRY (GrB_Vector_setElement_INT64 (*p_result, p [0], 0)) ;
            GRB_TRY (GrB_Vector_setElement_INT64 (*p_result, p [1], 1)) ;
        }
    }
    else
    {
        for (int64_t k = 0 ; k < n ; k++)
        {
            if (p [k] != -1)
            {
                GRB_TRY (GrB_Vector_setElement_FP64  (*x_result, v [k], k)) ;
                GRB_TRY (GrB_Vector_setElement_INT64 (*p_result, p [k], k)) ;
            }
        }
    }

    GRB_TRY (GrB_wait (*x_result, GrB_MATERIALIZE)) ;
    GRB_TRY (GrB_wait (*p_result, GrB_MATERIALIZE)) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

