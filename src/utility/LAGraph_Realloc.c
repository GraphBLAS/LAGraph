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

// The actual size_allocated on input can differ from nitems_old*size_of_item,
// and the size_allocated on output can be larger than nitems_new*size_of_item,
// as determined by the underlying memory manager.

// Usage:

//      p = LAGraph_Realloc (nitems_new, nitems_old, size_of_item, p,
//          &size_allocated, &ok)
//      if (ok)
//      {
//          p points to a block of at least nitems_new*size_of_item bytes and
//          the first part, of size min(nitems_new,nitems_old)*size_of_item,
//          has the same content as the old memory block if it was present.
//      }
//      else
//      {
//          p points to the old block, and size_allocated is left
//          unchanged.  This case never occurs if nitems_new < nitems_old.
//      }
//      on output, size_allocated is set to the actual size of the block of
//      memory

#include "LG_internal.h"

void *LAGraph_Realloc       // returns pointer to reallocated block of memory,
(                           // or original block if reallocation fails.
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // size of each item
    // input/output
    void *p,                // old block to reallocate
    size_t *size_allocated, // # of bytes actually allocated
    // output
    bool *ok                // true if successful, false otherwise
)
{

    //--------------------------------------------------------------------------
    // malloc a new block if p is NULL on input
    //--------------------------------------------------------------------------

    if (p == NULL)
    { 
        p = LAGraph_Malloc (nitems_new, size_of_item, size_allocated) ;
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

    if (!(*ok) || nitems_new > GxB_INDEX_MAX || size_of_item > GxB_INDEX_MAX)
    { 
        // overflow
        (*ok) = false ;
        return (NULL) ;
    }

    //--------------------------------------------------------------------------
    // reallocate an existing block to accommodate the change in # of items
    //--------------------------------------------------------------------------

    int64_t oldsize_allocated = (*size_allocated) ;

    //--------------------------------------------------------------------------
    // check for quick return
    //--------------------------------------------------------------------------

    if ((newsize == oldsize)
        || (newsize < oldsize && newsize >= oldsize_allocated/2)
        || (newsize > oldsize && newsize <= oldsize_allocated))
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
    size_t newsize_allocated = newsize ;

    if (LAGraph_Realloc_function == NULL)
    {

        //----------------------------------------------------------------------
        // use malloc/memcpy/free
        //----------------------------------------------------------------------

        // allocate the new space
        pnew = LAGraph_Malloc (nitems_new, size_of_item, &newsize_allocated) ;
        // copy over the data from the old block to the new block
        if (pnew != NULL)
        {
            // TODO: use a parallel memcpy
            memcpy (pnew, p, LAGraph_MIN (oldsize, newsize)) ;
            // free the old space
            LAGraph_Free (&p, oldsize_allocated) ;
        }
    }
    else
    {

        //----------------------------------------------------------------------
        // use realloc
        //----------------------------------------------------------------------

        pnew = LAGraph_Realloc_function (p, newsize_allocated) ;
    }

    //--------------------------------------------------------------------------
    // check if successful and return result
    //--------------------------------------------------------------------------

    if (pnew == NULL)
    {
        // realloc failed
        if (newsize < oldsize)
        { 
            // the attempt to reduce the size of the block failed, but the old
            // block is unchanged.  So pretend to succeed, but do not change
            // size_allocated since it must reflect the actual size of the
            // block.
            (*ok) = true ;
        }
        else
        { 
            // out of memory.  the old block is unchanged
            (*ok) = false ;
        }
    }
    else
    { 
        // realloc succeeded
        p = pnew ;
        (*ok) = true ;
        (*size_allocated) = newsize_allocated ;
    }

    return (p) ;
}

