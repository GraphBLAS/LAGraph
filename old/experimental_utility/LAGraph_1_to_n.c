//------------------------------------------------------------------------------
// LAGraph_1_to_n.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//****************************************************************************
// FIXME: this is not yet included in the test coverage suite

// Create either a GrB_INT64 or GrB_INT32 "ramp" vector 1:n

#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&v) ;         \
    LAGraph_Free ((void **)&I) ;                 \
    LAGraph_Free ((void **)&X) ;                 \
}

//****************************************************************************
/// @todo If this method gets promoted it should return GrB_Type for scalar
///       that is stored in the output vector.
///
GrB_Info LAGraph_1_to_n     // create an integer vector v = 1:n
(
    GrB_Vector *v_handle,   // vector to create
    GrB_Index n             // size of vector to create
)
{

    GrB_Info info ;
    GrB_Vector v = NULL ;
    int nthreads;
    LAGraph_GetNumThreads (&nthreads, NULL) ;
    nthreads = LAGraph_MIN (n / 4096, nthreads) ;
    nthreads = LAGraph_MAX (nthreads, 1) ;

    // allocate workspace
    GrB_Index *I = LAGraph_Malloc (n, sizeof (GrB_Index)) ;

    // create a 32-bit or 64-bit integer vector 1:n
    if (n > INT32_MAX)
    {
        int64_t *X = LAGraph_Malloc (n, sizeof (int64_t)) ;
        if (I == NULL || X == NULL)
        {
            LAGraph_FREE_ALL ;
            return (GrB_OUT_OF_MEMORY) ;
        }
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t k = 0 ; k < n ; k++)
        {
            I [k] = k ;
            X [k] = k+1 ;
        }
        LAGRAPH_OK (GrB_Vector_new (&v, GrB_INT64, n)) ;
        LAGRAPH_OK (GrB_Vector_build (v, I, X, n, GrB_PLUS_INT64)) ;
        LAGraph_Free ((void **)&X) ;  X = NULL;
    }
    else
    {
        int32_t *X = LAGraph_Malloc (n, sizeof (int32_t)) ;
        if (I == NULL || X == NULL)
        {
            LAGraph_FREE_ALL ;
            return (GrB_OUT_OF_MEMORY) ;
        }
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int64_t k = 0 ; k < n ; k++)
        {
            I [k] = k ;
            X [k] = k+1 ;
        }
        LAGRAPH_OK (GrB_Vector_new (&v, GrB_INT32, n)) ;
        LAGRAPH_OK (GrB_Vector_build (v, I, X, n, GrB_PLUS_INT32)) ;
        LAGraph_Free ((void **)&X) ;  X = NULL;
    }
    LAGraph_Free ((void **)&I) ;  I = NULL;

    // return result
    (*v_handle) = v ;  v = NULL;
    return (GrB_SUCCESS) ;
}
