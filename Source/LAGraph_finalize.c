//------------------------------------------------------------------------------
// LAGraph_finalize:  start LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

GrB_Info LAGraph_finalize ( )
{
    // free the complex type for LAGraph
    GrB_free (&LAGraph_Complex) ;

    // finalize GraphBLAS
    return (GrB_finalize ( )) :
}

