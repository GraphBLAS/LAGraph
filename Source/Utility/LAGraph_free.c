//------------------------------------------------------------------------------
// LAGraph_free:  wrapper for free
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Wrapper for free.

#include "LAGraph_internal.h"

void LAGraph_free
(
    void *p
)
{

    if (p != NULL)
    {
        LAGraph_free_function (p) ;
    }
}

