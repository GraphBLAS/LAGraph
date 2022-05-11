//------------------------------------------------------------------------------
// LAGraph_GetNumThreads: get the # of threads to use
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_get_nthreads: get # of threads that will be used by LAGraph.

#include "LG_internal.h"

int LAGraph_GetNumThreads
(
    // output:
    int *nthreads_hi,
    int *nthreads_lo,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (nthreads_hi != NULL && nthreads_lo != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // get number of threads
    //--------------------------------------------------------------------------

    (*nthreads_hi) = LG_nthreads_hi ;
    (*nthreads_lo) = LG_nthreads_lo ;
    return (GrB_SUCCESS) ;
}

