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
    GrB_free (&M) ;                 \
    GrB_free (&thunk) ;             \
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
    GrB_Matrix T = NULL, M = NULL, A = NULL ;
    GrB_Vector y = NULL, u = NULL, w = NULL ;
    GxB_Scalar thunk = NULL ;

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

GxB_set (GxB_BURBLE, true) ;    // FIXME: remove
    //--------------------------------------------------------------------------
    // count triangles: T<A> = A*A' using the plus_pair semiring
    //--------------------------------------------------------------------------

    // TODO: try computing just L = tril(T), which would be faster.
    // to do this, use M = tril(A) as the structural mask.
    // Then the computations below would be:
    //      y = sum (L) + sum (L')
    //      k = 2 * sum (L)
    //      centrality = (3*A*y - 2*(L*y + L'*y) + y) / k

    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_new (&T, GrB_FP64, n, n)) ;

#if 1
    GrB_TRY (GrB_Matrix_new (&M, GrB_FP64, n, n)) ;

    // M = tril (A,-1)
    GrB_TRY (GxB_Scalar_new (&thunk, GrB_INT64)) ;
    GrB_TRY (GxB_Scalar_setElement (thunk, (int64_t) (-1))) ;
    GrB_TRY (GxB_select (M, NULL, NULL, GxB_TRIL, A, thunk, NULL)) ;
    GrB_TRY (GrB_free (&thunk)) ;

    // T<M>= A*A'
    GrB_TRY (GrB_mxm (T, M, NULL, GxB_PLUS_PAIR_FP64, A, A, GrB_DESC_T1)) ;
    GrB_TRY (GrB_free (&M)) ;

    // y = sum (T'), where y(j) = sum (T (:,j)) and y(j)=0 if T(:,j) is empty
    GrB_TRY (GrB_Vector_new (&y, GrB_FP64, n)) ;
    GrB_TRY (GrB_assign (y, NULL, NULL, ((double) 0), GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_reduce (y, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, T,
        GrB_DESC_T0)) ;
    // y += sum (T)
    GrB_TRY (GrB_reduce (y, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, T,
        NULL)) ;

    // k = sum (y)
    double k = 0 ;
    GrB_TRY (GrB_reduce (&k, NULL, GrB_PLUS_MONOID_FP64, y, NULL)) ;

    // centrality = (3*A*y - 2* (T*y + T'*y) + y) / k

    // w = T*y
    GrB_TRY (GrB_Vector_new (&w, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (w, NULL, NULL, GxB_PLUS_SECOND_FP64, T, y, NULL)) ;
    // w += T'*y
    GrB_TRY (GrB_mxv (w, NULL, GrB_PLUS_FP64, GxB_PLUS_SECOND_FP64, T, y,
        GrB_DESC_T0)) ;

#else

    GrB_TRY (GrB_mxm (T, A, NULL, GxB_PLUS_PAIR_FP64, A, A, GrB_DESC_T1)) ;

    //--------------------------------------------------------------------------
    // y = sum (T), where y(i) = sum (T (i,:)) and y(i)=0 of T(i,:) is empty
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Vector_new (&y, GrB_FP64, n)) ;
    GrB_TRY (GrB_assign (y, NULL, NULL, ((double) 0), GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_reduce (y, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, T,
        NULL)) ;

    //--------------------------------------------------------------------------
    // k = sum (y)
    //--------------------------------------------------------------------------

    double k = 0 ;
    GrB_TRY (GrB_reduce (&k, NULL, GrB_PLUS_MONOID_FP64, y, NULL)) ;

    //--------------------------------------------------------------------------
    // centrality = (3*A*y - 2*T*y + y) / k
    //--------------------------------------------------------------------------

    // w = T*y
    GrB_TRY (GrB_Vector_new (&w, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (w, NULL, NULL, GxB_PLUS_SECOND_FP64, T, y, NULL)) ;

#endif

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

    // centrality = centrality / k
    GrB_TRY (GrB_apply (*centrality, NULL, NULL, GrB_TIMES_FP64,
        ((k == 0) ? 1.0 : (1.0/k)), *centrality, NULL)) ;

    // TODO: # of triangles is k/6, which could be returned as well
    printf ("# of triangles: %g\n", k/6) ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_WORK ;
GxB_set (GxB_BURBLE, false) ;   // FIXME: remove
    return (0) ;
}

