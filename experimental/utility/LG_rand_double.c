//------------------------------------------------------------------------------
// LAGraph_rand_double: return a random double
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

// LAGraph_rand_double: return a random double between 0 and 1, inclusive.
// Contributed by Tim Davis, Texas A&M

#include "LAGraph_internal.h"

double LAGraph_rand_double (uint64_t *seed)
{
    return (((double) LAGraph_rand64 (seed)) / ((double) UINT64_MAX)) ;
}
