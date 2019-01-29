//------------------------------------------------------------------------------
// LAGraph_init:  start LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

// a global value for returning the complex type in a Matrix Market file:
GrB_Type LAGraph_Complex = NULL ;

GrB_Info LAGraph_init ( )
{

    // initialize GraphBLAS
    GrB_Info info = GrB_init ( ) ;
    if (info != GrB_SUCCESS)
    {
        return (info) ;
    }

    // create the complex type for LAGraph
    return (GrB_Type_new (&LAGraph_Complex, sizeof (double complex))) ;
}

