//------------------------------------------------------------------------------
// LAGraph_randx: return a random double
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_randx: return a random double between 0 and 1, inclusive.

#include "LAGraph_internal.h"

double LAGraph_randx (uint64_t *seed)
{
    return (((double) LAGraph_rand (seed)) / ((double) UINT64_MAX)) ;
}

