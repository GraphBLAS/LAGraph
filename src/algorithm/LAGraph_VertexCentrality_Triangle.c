//------------------------------------------------------------------------------
// LAGraph_VertexCentrality_triangle: vertex triangle-centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University, ...

//------------------------------------------------------------------------------

// LAGraph_VertexCentrality_Triangle: computes the TriangleCentrality of
// an undirected graph.  No self edges are allowed on the input graph.

// P. Burkhardt, "Triangle centrality," https://arxiv.org/pdf/2105.00110.pdf,
// April 2021.

// TODO: this uses GxB* methods

//------------------------------------------------------------------------------

#include "LG_internal.h"

#define LAGRAPH_FREE_WORK           \
{                                   \
    GrB_free (&T) ;                 \
    GrB_free (&u) ;                 \
    GrB_free (&w) ;                 \
    GrB_free (&y) ;                 \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    GrB_free (centrality) ;         \
}

//------------------------------------------------------------------------------
// LAGraph_VertexCentrality_Triangle: vertex triangle-centrality
//------------------------------------------------------------------------------

int LAGraph_VertexCentrality_Triangle       // vertex triangle-centrality
(
    // outputs:
    GrB_Vector *centrality,     // centrality(i): triangle centrality of i
    // inputs:
    LAGraph_Graph G,            // input graph
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix T = NULL, A = NULL ;
    GrB_Vector y = NULL, u = NULL, w = NULL ;

    LG_CHECK (centrality == NULL, -1, "centrality is NULL") ;
    (*centrality) = NULL ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_pattern_is_symmetric == LAGRAPH_TRUE))
    {
        // the pattern of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -105, "G->A must be symmetric") ;
    }

    // no self edges can be present
    LG_CHECK (G->ndiag != 0, -104, "G->ndiag must be zero") ;

    //--------------------------------------------------------------------------
    // count triangles: T<A> = A*A' using the plus_pair semiring
    //--------------------------------------------------------------------------

    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_new (&T, GrB_FP64, n, n)) ;
    GrB_TRY (GrB_mxm (T, A, NULL, GxB_PLUS_PAIR_FP64, A, A, GrB_DESC_T1)) ;

    //--------------------------------------------------------------------------
    // y = sum (T), where y(i) = sum (T (i,:)) and y(i)=0 of T(i,:) is empty
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&y, GrB_FP64, n)) ;
    GrB_TRY (GrB_assign (y, NULL, NULL, ((double) 0), GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_reduce (y, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, T,
        NULL)) ;

    //--------------------------------------------------------------------------
    // ntriangles = sum (y), the total number of triangles
    //--------------------------------------------------------------------------

    double ntriangles = 0 ;
    GrB_TRY (GrB_reduce (&ntriangles, NULL, GrB_PLUS_MONOID_FP64, y, NULL)) ;

    //--------------------------------------------------------------------------
    // centrality = (3*A*y - 2*T*y + y) / ntriangles
    //--------------------------------------------------------------------------

    // w = T*y
    GrB_TRY (GrB_Vector_new (&w, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (w, NULL, NULL, GxB_PLUS_SECOND_FP64, T, y, NULL)) ;

    // w = (-2)*w
    double minus_two = -2 ;
    GrB_TRY (GrB_apply (w, NULL, NULL, GrB_TIMES_FP64, minus_two, w, NULL)) ;

    // u = A*y
    GrB_TRY (GrB_Vector_new (&u, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (u, NULL, NULL, GxB_PLUS_SECOND_FP64, A, y, NULL)) ;

    // u = 3*u
    double three = 3 ;
    GrB_TRY (GrB_apply (u, NULL, NULL, GrB_TIMES_FP64, three, u, NULL)) ;

    // centrality = (u + w + y)
    GrB_TRY (GrB_Vector_dup (centrality, y)) ;
    GrB_TRY (GrB_eWiseAdd (*centrality, NULL, GrB_PLUS_FP64, GrB_PLUS_FP64,
        u, w, NULL)) ;

    // centrality = centrality / ntriangles
    GrB_TRY (GrB_apply (*centrality, NULL, NULL, GrB_TIMES_FP64,
        ((ntriangles == 0) ? 1.0 : (1.0/ntriangles)), *centrality, NULL)) ;

    // TODO: also return # of triangles

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_WORK ;
    return (0) ;
}

