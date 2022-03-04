//------------------------------------------------------------------------------
// LAGraph_KindName: return the name of a kind
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_KindName
(
    char *name,     // name of the kind (user provided array of size at least
                    // LAGRAPH_MAX_NAME_LEN)
    LAGraph_Kind kind,  // graph kind
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
        case LAGraph_ADJACENCY_UNDIRECTED : strcpy (name, "undirected"); break ;
        case LAGraph_ADJACENCY_DIRECTED :   strcpy (name, "directed")  ; break ;
        case LAGraph_KIND_UNKNOWN :         strcpy (name, "unknown")   ; break ;
        default : LG_ASSERT_MSG (false, GrB_INVALID_VALUE, "invalid kind") ;
    }

    return (GrB_SUCCESS) ;
}

