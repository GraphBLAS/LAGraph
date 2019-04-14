//------------------------------------------------------------------------------
// LAGraph_malloc:  wrapper for malloc
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Wrapper for malloc.

// TODO also need a wrapper for calloc and realloc.

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// global space
//------------------------------------------------------------------------------

void * (* LAGraph_malloc_function  ) (size_t)         = malloc ;
void * (* LAGraph_calloc_function  ) (size_t, size_t) = calloc ;
void * (* LAGraph_realloc_function ) (void *, size_t) = realloc ;
void   (* LAGraph_free_function    ) (void *)         = free ;

bool LAGraph_malloc_is_thread_safe =
    #ifdef MATLAB_MEX_FILE
        false       // mxMalloc is not thread-safe
    #else
        true        // ANCI C malloc, TBB scalable_malloc, etc are thread safe
    #endif
        ;

//------------------------------------------------------------------------------
// LAGraph_malloc
//------------------------------------------------------------------------------

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
    return (LAGraph_malloc_function (nitems * size_of_item)) ;
}

