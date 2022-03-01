//------------------------------------------------------------------------------
// LAGraph_SetNumThreads: set the # of threads to use
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_SetNumThreads
(
    // input:
    int nthreads,       // # of threads to use
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    //--------------------------------------------------------------------------
    // set number of threads to use
    //--------------------------------------------------------------------------

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

    return (GrB_SUCCESS) ;
}

