//------------------------------------------------------------------------------
// LAGraph_KindName: return the name of a kind
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_KindName        // returns 0 if successful, -1 if failure
(
    char *name,             // name of the kind
    LAGraph_Kind kind,      // graph kind
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the name of the kind
    //--------------------------------------------------------------------------

    switch (kind)
    {
        case LAGRAPH_ADJACENCY_UNDIRECTED : strcpy (name, "undirected"); break ;
        case LAGRAPH_ADJACENCY_DIRECTED :   strcpy (name, "directed")  ; break ;
        case LAGRAPH_KIND_UNKNOWN :         strcpy (name, "unknown")   ; break ;
        default : 
            LG_ASSERT_MSG (false, GrB_INVALID_VALUE, "invalid kind") ; // RETVAL
    }

    return (GrB_SUCCESS) ;
}

