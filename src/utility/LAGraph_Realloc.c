//------------------------------------------------------------------------------
// LAGraph_Realloc: wrapper for realloc
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// If p is non-NULL on input, it points to a previously allocated object of
// size at least nitems_old * size_of_item.  The object is reallocated to be of
// size at least nitems_new * size_of_item.  If p is NULL on input, then a new
// object of that size is allocated.  On success, a pointer to the new object
// is returned, and ok is returned as true.  If the allocation fails, ok is set
// to false and a pointer to the old (unmodified) object is returned.

// Usage:

//      p = LAGraph_Realloc (nitems_new, nitems_old, size_of_item, p, &ok)
//      if (ok)
//      {
//          p points to a block of at least nitems_new*size_of_item bytes and
//          the first part, of size min(nitems_new,nitems_old)*size_of_item,
//          has the same content as the old memory block if it was present.
//      }
//      else
//      {
//          p points to the old block, unchanged.  This case never occurs if
//          nitems_new < nitems_old.
//      }

#include "LG_internal.h"

void *LAGraph_Realloc       // returns pointer to reallocated block of memory,
(                           // or original block if reallocation fails.
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // size of each item
    // input/output
    void *p,                // old block to reallocate
    // output
    bool *ok                // true if successful, false otherwise
)
{

    //--------------------------------------------------------------------------
    // malloc a new block if p is NULL on input
    //--------------------------------------------------------------------------

    if (p == NULL)
    {
        p = LAGraph_Malloc (nitems_new, size_of_item) ;
        (*ok) = (p != NULL) ;
        return (p) ;
    }

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // make sure at least one item is allocated
    nitems_old = LAGraph_MAX (1, nitems_old) ;
    nitems_new = LAGraph_MAX (1, nitems_new) ;

    // make sure at least one byte is allocated
    size_of_item = LAGraph_MAX (1, size_of_item) ;

    size_t newsize, oldsize ;
    (*ok) = LG_Multiply_size_t (&newsize, nitems_new, size_of_item)
         && LG_Multiply_size_t (&oldsize, nitems_old, size_of_item) ;

    if (!(*ok) || nitems_new > GrB_INDEX_MAX || size_of_item > GrB_INDEX_MAX)
    {
        // overflow
        (*ok) = false ;
        return (NULL) ;
    }

    //--------------------------------------------------------------------------
    // reallocate an existing block to accommodate the change in # of items
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // check for quick return
    //--------------------------------------------------------------------------

    if (newsize == oldsize)
    {
        // If the block does not change, or is shrinking but only by a small
        // amount, or is growing but still fits inside the existing block,
        // then leave the block as-is.
        (*ok) = true ;
        return (p) ;
    }

    //--------------------------------------------------------------------------
    // reallocate the memory, or use malloc/memcpy/free
    //--------------------------------------------------------------------------

    void *pnew = NULL ;

    if (LAGraph_Realloc_function == NULL)
    {

        //----------------------------------------------------------------------
        // use malloc/memcpy/free
        //----------------------------------------------------------------------

        // allocate the new space
        pnew = LAGraph_Malloc (nitems_new, size_of_item) ;
        // copy over the data from the old block to the new block
        if (pnew != NULL)
        {
            // use a parallel memcpy
            memcpy (pnew, p, LAGraph_MIN (oldsize, newsize)) ;
            // free the old space
            LAGraph_Free (&p) ;
        }
    }
    else
    {

        //----------------------------------------------------------------------
        // use realloc
        //----------------------------------------------------------------------

        pnew = LAGraph_Realloc_function (p, newsize) ;
    }

    //--------------------------------------------------------------------------
    // check if successful and return result
    //--------------------------------------------------------------------------

    // If the attempt to reduce the size of the block failed, the old block is
    // unchanged.  So pretend to succeed.

    (*ok) = (newsize < oldsize) || (pnew != NULL) ;
    p = (pnew != NULL) ? pnew : p ;
    return (p) ;
}
