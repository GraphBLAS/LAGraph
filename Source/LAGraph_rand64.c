//------------------------------------------------------------------------------
// LAGraph_rand64: return a random 64-bit integer
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_rand64: return a random 64-bit unsigned integer.

#include "LAGraph_internal.h"

uint64_t LAGraph_rand64 (uint64_t *seed)
{
    uint64_t i = 0 ;
    for (int k = 0 ; k < 5 ; k++)
    {
        i = LAGRAPH_RAND_MAX * i + LAGraph_rand (seed) ;
    }
    return (i) ;
}

