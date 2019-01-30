//------------------------------------------------------------------------------
// LAGraph_free:  wrapper for free
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

void LAGraph_free
(
    void **p                // *p is freed and set to NULL
)
{

    if (p != NULL)
    {
        if (*p != NULL)
        {
            free (p) ;
            (*p) = NULL ;
        }
    }
}

