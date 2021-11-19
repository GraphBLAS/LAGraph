//------------------------------------------------------------------------------
// LAGraph_rand64: return a random 64-bit integer
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite
// FIXME: remove this and use LAGraph_Random.

// LAGraph_rand64: return a random 64-bit unsigned integer.
// Contributed by Tim Davis, Texas A&M.

#include <stdint.h>
#include "LG_internal.h"

#define LAGRAPH_RAND_MAX 32767

uint64_t LAGraph_rand64 (uint64_t *seed)
{
    uint64_t i = 0 ;
    for (int k = 0 ; k < 5 ; k++)
    {
        i = LAGRAPH_RAND_MAX * i + LAGraph_Random15 (seed) ;
    }
    return (i) ;
}
