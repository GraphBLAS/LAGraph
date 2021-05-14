//------------------------------------------------------------------------------
// LAGraph_SetNumThreads: set the # of threads to use
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_SetNumThreads       // returns 0 if successful, -1 if failure
(
    int nthreads,
    char *msg
)
{

    LG_CLEAR_MSG ;

    #if LG_SUITESPARSE
    {
        // SuiteSparse:GraphBLAS: set # of threads with global setting
        GrB_TRY (GxB_set (GxB_NTHREADS, nthreads)) ;
    }
    #elif defined ( _OPENMP )
    {
        // set # of threads with OpenMP global setting
        omp_set_num_threads (nthreads) ;
    }
    #else
    {
        // nothing to do
    }
    #endif

    return (0) ;
}

