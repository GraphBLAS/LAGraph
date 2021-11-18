//------------------------------------------------------------------------------
// LAGraph_Calloc:  wrapper for calloc
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

void *LAGraph_Calloc
(
    size_t nitems,          // number of items
    size_t size_of_item     // size of each item
)
{

    // make sure at least one item is allocated
    nitems = LAGraph_MAX (1, nitems) ;

    // make sure at least one byte is allocated
    size_of_item = LAGraph_MAX (1, size_of_item) ;

    // compute the size and check for integer overflow
    size_t size ;
    bool ok = LG_Multiply_size_t (&size, nitems, size_of_item) ;
    if (!ok || nitems > GrB_INDEX_MAX || size_of_item > GrB_INDEX_MAX)
    {
        // overflow
        return (NULL) ;
    }

    // calloc the space
    void *p = NULL ;
    if (LAGraph_Calloc_function != NULL)
    {
        // use the calloc function
        p = LAGraph_Calloc_function (nitems, size_of_item) ;
    }
    else
    {
        // calloc function not available; use malloc and memset
        p = LAGraph_Malloc (nitems, size_of_item) ;
        if (p != NULL)
        {
            memset (p, 0, size) ;
        }
    }
    return (p) ;
}
