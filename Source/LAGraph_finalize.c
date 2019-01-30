//------------------------------------------------------------------------------
// LAGraph_finalize:  finish LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Finalize LAGraph and then finilize GraphBLAS.

#include "LAGraph_internal.h"

GrB_Info LAGraph_finalize ( )
{
    // free the complex type and operators for LAGraph
    LAGraph_free_global ( ) ;

    // finalize GraphBLAS
    return (GrB_finalize ( )) ;
}

