//------------------------------------------------------------------------------
// LAGraph_SetNumThreads: set the # of threads to use
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_SetNumThreads
(
    // input:
    int nthreads_hi,
    int nthreads_lo,
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

    nthreads_hi = LAGRAPH_MAX (nthreads_hi, 1) ;
    nthreads_lo = LAGRAPH_MAX (nthreads_lo, 1) ;

    #if LAGRAPH_SUITESPARSE
    {
        // SuiteSparse:GraphBLAS: set # of threads with global setting
        GRB_TRY (GxB_set (GxB_NTHREADS, nthreads_lo)) ;
    }
    #endif

    // set # of LAGraph threads
    LG_nthreads_hi = nthreads_hi ;      // for LAGraph itself, if nested
                                        // regions call GraphBLAS
    LG_nthreads_lo = nthreads_lo ;      // for lower-level parallelism

    return (GrB_SUCCESS) ;
}

