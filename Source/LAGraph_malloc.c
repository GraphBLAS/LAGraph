//------------------------------------------------------------------------------
// LAGraph_malloc:  wrapper for malloc
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

GrB_Info LAGraph_malloc
(
    void **p,               // pointer to allocated block of memory
    size_t nitems,          // number of items
    size_t size_of_item     // size of each item
)
{
    if (p == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    // make sure at least one item is allocated
    nitems = LAGRAPH_MAX (1, nitems) ;

    // make sure at least one byte is allocated
    size_of_item = LAGRAPH_MAX (1, size_of_item) ;

    // TODO check for integer overflow
    (*p) = malloc (nitems * size_of_item) ;

    return (((*p) == NULL) ? GrB_OUT_OF_MEMORY : GrB_SUCCESS) ;
}

