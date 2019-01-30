//------------------------------------------------------------------------------
// LAGraph_rand: return a random number
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// LAGraph_rand: A simple thread-safe random number generator that returns a
// random number between 0 and LAGRAPH_RAND_MAX.  The quality of the random
// values it generates is very low, but this is not important.  This method is
// used to create random test matrices, which must be identical on any
// operating system.

#include "LAGraph_internal.h"

uint64_t LAGraph_rand (uint64_t *seed)
{
   (*seed) = (*seed) * 1103515245 + 12345 ;
   return (((*seed)/65536) % (LAGRAPH_RAND_MAX + 1)) ;
}

