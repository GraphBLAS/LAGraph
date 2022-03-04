//------------------------------------------------------------------------------
// LG_Random.c: simple and portable random number generator
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

// return a random number between 0 and LG_RANDOM15_MAX
GrB_Index LG_Random15 (uint64_t *seed)
{
   (*seed) = (*seed) * 1103515245 + 12345 ;
   return (((*seed) / 65536) % (LG_RANDOM15_MAX + 1)) ;
}

// return a random uint64_t, in range 0 to LG_RANDOM60_MAX
GrB_Index LG_Random60 (uint64_t *seed)
{
    GrB_Index i ;
    i = LG_Random15 (seed) ;
    i = LG_Random15 (seed) + LG_RANDOM15_MAX * i ;
    i = LG_Random15 (seed) + LG_RANDOM15_MAX * i ;
    i = LG_Random15 (seed) + LG_RANDOM15_MAX * i ;
    i = i % (LG_RANDOM60_MAX + 1) ;
    return (i) ;
}

