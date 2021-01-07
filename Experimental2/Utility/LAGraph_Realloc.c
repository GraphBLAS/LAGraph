//------------------------------------------------------------------------------
// LAGraph_Realloc: wrapper for realloc
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// If p is non-NULL on input, it points to a previously allocated object of
// size nitems_old * size_of_item.  The object is reallocated to be of size
// nitems_new * size_of_item.  If p is NULL on input, then a new object of that
// size is allocated.  On success, a pointer to the new object is returned, and
// ok is returned as true.  If the allocation fails, ok is set to false and a
// pointer to the old (unmodified) object is returned.

// Usage:

//      p = LAGraph_Realloc (nnew, nold, size, p, &ok)
//      if (ok)
//      {
//          p points to a space of size at least nnew*size, and the first
//          part, of size min(nnew,nold)*size, has the same content as
//          the old memory space if it was present.
//      }
//      else
//      {
//          p points to the old space of size nold*size, which is left
//          unchanged.  This case never occurs if nnew < nold.
//      }

#include "LAGraph_Internal.h"

void *LAGraph_Realloc       // returns pointer to reallocated block of memory,
(                           // or original block if reallocation fails.
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // size of each item
    void *p,                // old object to reallocate
    bool *ok                // true if successful, false otherwise
)
{
    size_t size, oldsize ;

    // make sure at least one item is allocated
    nitems_old = LAGraph_MAX (1, nitems_old) ;
    nitems_new = LAGraph_MAX (1, nitems_new) ;

    // make sure at least one byte is allocated
    size_of_item = LAGraph_MAX (1, size_of_item) ;

    (*ok) = LAGraph_Multiply_size_t (&size,    nitems_new, size_of_item)
         && LAGraph_Multiply_size_t (&oldsize, nitems_old, size_of_item) ;

    if (!(*ok) || nitems_new > GxB_INDEX_MAX || size_of_item > GxB_INDEX_MAX)
    { 
        // overflow
        (*ok) = false ;
    }
    else if (p == NULL)
    { 
        // a fresh object is being allocated
        p = (void *) LAGraph_Malloc (nitems_new, size_of_item) ;
        (*ok) = (p != NULL) ;
    }
    else if (nitems_old == nitems_new)
    { 
        // the object does not change; do nothing
        (*ok) = true ;
    }
    else
    {
        // change the size of the object from nitems_old to nitems_new
        void *pnew ;

        //----------------------------------------------------------------------
        // reallocate the memory
        //----------------------------------------------------------------------

        if (LAGraph_Realloc_function != NULL)
        {

            //------------------------------------------------------------------
            // use realloc
            //------------------------------------------------------------------

            pnew = (void *) LAGraph_Realloc_function (p, size) ;

        }
        else
        {

            //------------------------------------------------------------------
            // no realloc function: mimic with malloc and memcpy
            //------------------------------------------------------------------

            // malloc the new space
            pnew = (void *) LAGraph_Malloc_function (size) ;

            // copy over the data from the old space to the new space
            if (pnew != NULL)
            {
                // TODO: use a parallel memcpy
                memcpy (pnew, p, LAGraph_MIN (oldsize, size)) ;

                // free the old space
                LAGraph_Free_function (p) ;
            }
        }

        //----------------------------------------------------------------------
        // check if successful
        //----------------------------------------------------------------------

        if (pnew == NULL)
        {
            if (nitems_new < nitems_old)
            {
                // the attempt to reduce the size of the block failed, but
                // the old block is unchanged.  So pretend to succeed.
                (*ok) = true ;
            }
            else
            { 
                // out of memory
                (*ok) = false ;
            }
        }
        else
        { 
            // success
            p = pnew ;
            (*ok) = true ;
        }
    }

    //--------------------------------------------------------------------------
    // return newly allocated block, or old block if reallocation failed
    //--------------------------------------------------------------------------

    return (p) ;
}

