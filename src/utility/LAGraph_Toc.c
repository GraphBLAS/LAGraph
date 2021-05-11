//------------------------------------------------------------------------------
// LAGraph_Toc: return time since last call to LAGraph_Tic
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Toc             // 0 if successful, -1 on failure
(
    double *t,              // time since last call to LAGraph_Tic
    const double tic [2],   // tic from last call to LAGraph_Tic
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (t == NULL, -1, "&t is NULL") ;
    LG_CHECK (tic == NULL, -1, "tic is NULL") ;

    //--------------------------------------------------------------------------
    // get the time since the last call to LAGraph_Toc
    //--------------------------------------------------------------------------

    double toc [2] ;
    LAGraph_TRY (LAGraph_Tic (toc, msg)) ;
    (*t) = (toc [0] - tic [0]) + 1e-9 * (toc [1] - tic [1]) ;

    return (0) ;
}

