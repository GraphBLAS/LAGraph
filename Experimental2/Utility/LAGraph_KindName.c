//------------------------------------------------------------------------------
// LAGraph_KindName: return the name of a kind
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#include "LAGraph_Internal.h"

int LAGraph_KindName        // returns 0 if successful, -1 if failure
(
    char **name,            // name of the kind
    LAGraph_Kind kind,      // graph kind
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LAGraph_CLEAR_MSG ;
    LAGraph_CHECK (name == NULL, -1, "name is NULL") ;

    //--------------------------------------------------------------------------
    // determine the name of the kind
    //--------------------------------------------------------------------------

    switch (kind)
    {
        case LAGRAPH_ADJACENCY_UNDIRECTED : (*name) = "undirected" ; break ;
        case LAGRAPH_ADJACENCY_DIRECTED :   (*name) = "directed"   ; break ;
        case LAGRAPH_KIND_UNKNOWN :         (*name) = "unknown"    ; break ;
        default : LAGraph_CHECK (false, -1, "invalid kind") ;
    }

    return (0) ;
}

