//------------------------------------------------------------------------------
// LAGraph_GetNumThreads: get the # of threads to use
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_get_nthreads: get # of threads that will be used by LAGraph.

#include "LG_internal.h"

int LAGraph_GetNumThreads   // returns 0 if successful, or -1 if failure
(
    int *nthreads,          // # of threads to use
    char *msg
)
{

    LG_CLEAR_MSG ;
    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    {
        // SuiteSparse:GraphBLAS: get # of threads from global setting
        GrB_TRY (GxB_get (GxB_NTHREADS, nthreads)) ;
    }
    #elif defined ( _OPENMP )
    {
        // get # of theads from OpenMP global setting
        (*nthreads) = omp_get_max_threads ( ) ;
    }
    #else
    {
        // single-threade if not using SuiteSparse:GraphBLAS or OpenMP
        (*nthreads) = 1 ;
    }
    #endif
    return (0) ;
}

