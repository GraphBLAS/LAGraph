//------------------------------------------------------------------------------
// LAGraph_Toc: return time since last call to LAGraph_Tic
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Toc
(
    // output:
    double *t,              // time since last call to LAGraph_Tic
    // input:
    const double tic [2],   // tic from last call to LAGraph_Tic
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT_MSG (t != NULL, GrB_NULL_POINTER, "&t != NULL") ;
    LG_ASSERT (tic != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // get the time since the last call to LAGraph_Toc
    //--------------------------------------------------------------------------

    double toc [2] ;
    LG_TRY (LAGraph_Tic (toc, msg)) ;
    (*t) = (toc [0] - tic [0]) + 1e-9 * (toc [1] - tic [1]) ;

    return (GrB_SUCCESS) ;
}

