//------------------------------------------------------------------------------
// LAGraph_Calloc:  wrapper for calloc
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_Internal.h"

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
    bool ok = LAGraph_Multiply_size_t (&size, nitems, size_of_item) ;
    if (!ok || nitems > GxB_INDEX_MAX || size_of_item > GxB_INDEX_MAX)
    {
        // overflow
        return (NULL) ;
    }

    // calloc the space
    if (LAGraph_Calloc_function != NULL)
    {
        // use the calloc function
        return (LAGraph_Calloc_function (nitems, size_of_item)) ;
    }
    else
    {
        // calloc function not available; use malloc and memset
        void *p = LAGraph_Malloc_function (size) ;
        if (p != NULL)
        {
            memset (p, 0, size) ;
        }
        return (p) ;
    }
}

