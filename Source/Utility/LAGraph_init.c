//------------------------------------------------------------------------------
// LAGraph_init:  start LAGraph
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Initialize GraphBLAS, and then LAGraph.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL                    \
{                                           \
    LAGraph_free_global ( ) ;               \
}

GrB_Info LAGraph_init ( )
{
    // initialize GraphBLAS
    LAGRAPH_OK (GrB_init (GrB_NONBLOCKING)) ;

    // allocate all global objects
    LAGRAPH_OK (LAGraph_alloc_global ( )) ;

    return (GrB_SUCCESS) ;
}

