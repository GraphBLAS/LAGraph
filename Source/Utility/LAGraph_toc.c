//------------------------------------------------------------------------------
// LAGraph_toc:  a portable timer for accurate performance measurements
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// LAGraph_toc: return the time since the last LAGraph_tic
//------------------------------------------------------------------------------

double LAGraph_toc          // returns time since last LAGraph_tic
(
    const double tic [2]    // tic from last call to LAGraph_tic
)
{

    double toc [2] ;
    LAGraph_tic (toc) ;
    return ((toc [0] - tic [0]) + 1e-9 * (toc [1] - tic [1])) ;
}

