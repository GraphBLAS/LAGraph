//------------------------------------------------------------------------------
// LAGraph_New:  create a new graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// If succesful, the matrix A is "moved" into G->A, and the caller's A is set
// to NULL.

#include "LG_internal.h"

int LAGraph_New         // returns 0 if successful, -1 if failure
(
    LAGraph_Graph *G,   // the graph to create, NULL if failure
    GrB_Matrix *A,      // the adjacency matrix of the graph, may be NULL
    LAGraph_Kind kind,  // the kind of graph, may be LAGRAPH_UNKNOWN
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_CHECK (G == NULL, -1, "&G cannot be NULL on input") ;

    //--------------------------------------------------------------------------
    // allocate the graph
    //--------------------------------------------------------------------------

    (*G) = LAGraph_Malloc (1, sizeof (struct LAGraph_Graph_struct)) ;
    LG_CHECK (*G == NULL, -1, "out of memory") ;

    //--------------------------------------------------------------------------
    // assign its primary components
    //--------------------------------------------------------------------------

    if (A == NULL || *A == NULL)
    {
        // G is incomplete, to be finalized later
        (*G)->A = NULL ;
    }
    else
    {
        // move &A into the graph and set &A to NULL to denote to the caller
        // that it is now a component of G.  The graph G is not opaque, so the
        // caller can get A back with A = G->A, but this helps with memory
        // management, since LAGraph_Delete (&G,msg) frees G->A, and if the
        // caller also does GrB_free (&A), a double-free would occur if this
        // move does not set A to NULL.
        (*G)->A = (*A) ;
        (*A) = NULL ;
    }
    (*G)->kind = kind ;

    //--------------------------------------------------------------------------
    // clear its cached properties
    //--------------------------------------------------------------------------

    (*G)->AT = NULL ;
    (*G)->rowdegree = NULL ;
    (*G)->coldegree = NULL ;
    (*G)->A_pattern_is_symmetric =
        (kind == LAGRAPH_ADJACENCY_UNDIRECTED) ? LAGRAPH_TRUE : LAGRAPH_UNKNOWN ;
    (*G)->ndiag = LAGRAPH_UNKNOWN ;

    return (0) ;
}

