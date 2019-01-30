//------------------------------------------------------------------------------
// LAGraph_free:  wrapper for free
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Wrapper for free.

// TODO: should be able to modify the malloc/free routines used; also in
// GraphBLAS.

#include "LAGraph_internal.h"

void LAGraph_free
(
    void *p
)
{

    if (p != NULL)
    {
        free (p) ;
    }
}

