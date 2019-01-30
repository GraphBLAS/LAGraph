//------------------------------------------------------------------------------
// LAGraph_init:  start LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL                    \
{                                           \
    LAGraph_free_global ( ) ;               \
}

GrB_Info LAGraph_init ( )
{
    // initialize GraphBLAS
    LAGRAPH_OK (GrB_init ( )) ;

    // allocate all global objects
    LAGRAPH_OK (LAGraph_alloc_global ( )) ;
}

