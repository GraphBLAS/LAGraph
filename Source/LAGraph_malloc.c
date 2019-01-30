//------------------------------------------------------------------------------
// LAGraph_malloc:  wrapper for malloc
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Wrapper for malloc.

#include "LAGraph_internal.h"

void *LAGraph_malloc
(
    size_t nitems,          // number of items
    size_t size_of_item     // size of each item
)
{

    // make sure at least one item is allocated
    nitems = LAGRAPH_MAX (1, nitems) ;

    // make sure at least one byte is allocated
    size_of_item = LAGRAPH_MAX (1, size_of_item) ;

    // check for integer overflow
    if ((double) nitems * (double) size_of_item > (double) INT64_MAX)
    {
        return (NULL) ;
    }

    // malloc the space
    return (malloc (nitems * size_of_item)) ;
}

