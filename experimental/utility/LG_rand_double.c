//------------------------------------------------------------------------------
// LAGraph_rand_double: return a random double
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite
// FIXME: remove this and use LAGraph_Random.

// LAGraph_rand_double: return a random double between 0 and 1, inclusive.
// Contributed by Tim Davis, Texas A&M

#include <limits.h>
#include <LAGraphX.h>

double LAGraph_rand_double (uint64_t *seed)
{
    return (((double) LAGraph_rand64 (seed)) / ((double) UINT64_MAX)) ;
}
