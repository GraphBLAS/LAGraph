//------------------------------------------------------------------------------
// LG_Random.c: simple and portable random number generator
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LG_internal.h"

// return a random number between 0 and LAGRAPH_RANDOM15_MAX
GrB_Index LAGraph_Random15 (uint64_t *seed)
{
   (*seed) = (*seed) * 1103515245 + 12345 ;
   return (((*seed) / 65536) % (LAGRAPH_RANDOM15_MAX + 1)) ;
}

// return a random uint64_t, in range 0 to LAGRAPH_RANDOM60_MAX
GrB_Index LAGraph_Random60 (uint64_t *seed)
{
    GrB_Index i ;
    i = LAGraph_Random15 (seed) ;
    i = LAGraph_Random15 (seed) + LAGRAPH_RANDOM15_MAX * i ;
    i = LAGraph_Random15 (seed) + LAGRAPH_RANDOM15_MAX * i ;
    i = LAGraph_Random15 (seed) + LAGRAPH_RANDOM15_MAX * i ;
    i = i % (LAGRAPH_RANDOM60_MAX + 1) ;
    return (i) ;
}

