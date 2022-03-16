//------------------------------------------------------------------------------
// LAGraph_New:  create a new graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// If succesful, the matrix A is moved into G->A, and the caller's A is set
// to NULL.

#include "LG_internal.h"

int LAGraph_New
(
    // output:
    LAGraph_Graph *G,   // the graph to create, NULL if failure
    // input/output:
    GrB_Matrix    *A,   // the adjacency matrix of the graph, may be NULL.
                        // A is moved into G as G->A, and A itself is set
                        // to NULL to denote that is now a part of G.
                        // That is, { G->A = A ; A = NULL ; } is performed.
                        // When G is deleted, G->A is freed.  If A is NULL,
                        // the graph is invalid until G->A is set.
    // input:
    LAGraph_Kind kind,  // the kind of graph. This may be LAGRAPH_UNKNOWN,
                        // which must then be revised later before the
                        // graph can be used.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (G != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // allocate the graph
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Malloc ((void **) G, 1, sizeof (struct LAGraph_Graph_struct), msg)) ;

    //--------------------------------------------------------------------------
    // initialize its members
    //--------------------------------------------------------------------------

    (*G)->A = NULL ;
    (*G)->kind = LAGraph_KIND_UNKNOWN ;
    (*G)->AT = NULL ;
    (*G)->rowdegree = NULL ;
    (*G)->coldegree = NULL ;
    (*G)->structure_is_symmetric = LAGRAPH_UNKNOWN ;
    (*G)->ndiag = LAGRAPH_UNKNOWN ;
    (*G)->emin = NULL ;
    (*G)->emin_kind = LAGRAPH_UNKNOWN ;
    (*G)->emax = NULL ;
    (*G)->emax_kind = LAGRAPH_UNKNOWN ;
//  (*G)->nonzero = LAGRAPH_UNKNOWN ;       // future property

    //--------------------------------------------------------------------------
    // assign its primary components
    //--------------------------------------------------------------------------

    if ((A != NULL) && (*A != NULL))
    {
        // move &A into the graph and set &A to NULL to denote to the caller
        // that it is now a component of G.  The graph G is not opaque, so the
        // caller can get A back with A = G->A, but this helps with memory
        // management, since LAGraph_Delete (&G,msg) frees G->A, and if the
        // caller also does GrB_free (&A), a double-free would occur if this
        // move does not set A to NULL.
        (*G)->A = (*A) ;
        (*A) = NULL ;

        (*G)->kind = kind ;
        (*G)->structure_is_symmetric =
            (kind == LAGraph_ADJACENCY_UNDIRECTED)
            ? LAGraph_TRUE
            : LAGRAPH_UNKNOWN ;
    }

    return (GrB_SUCCESS) ;
}
