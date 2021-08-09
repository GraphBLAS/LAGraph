//------------------------------------------------------------------------------
// LAGraph_VertexCentrality_triangle: vertex triangle-centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University, ...

//------------------------------------------------------------------------------

// LAGraph_VertexCentrality_Triangle: computes the TriangleCentrality of
// an undirected graph.  No self edges are allowed on the input graph.
// Methods 2 and 3 can tolerate any edge weights (they are ignored; only the
// pattern of G->A is used).  Methods 1 and 1.5 require unit edge weights
// (this could be modified); results are undefined if this condition doesn't
// hold.

// P. Burkhardt, "Triangle centrality," https://arxiv.org/pdf/2105.00110.pdf,
// April 2021.

// Methods 2 and 3 require SuiteSparse:GraphBLAS.  Method 3 is by far the
// fastest.

// TC1: in python
//
//      def triangle_centrality1(A):
//          T = A.mxm(A, mask=A)
//          y = T.reduce_vector()
//          k = y.reduce_float()
//          return(1/k)*(3*(A @ y) - 2*(T @ y) + y)
//          note: T@y is wrong. should be plus_second semiring

// TC1.5: in python
//
//      def triangle_centrality1(A):
//          T = A.mxm(A, mask=A, desc=descriptor.ST1)
//          y = T.reduce_vector()
//          k = y.reduce_float()
//          return(1/k)*(3*(A @ y) - 2*(T @ y) + y)
//          note: T@y is wrong. should be plus_second semiring

// TC2: in python
//
//      def triangle_centrality2(A):
//          T = A.plus_pair(A, mask=A, desc=descriptor.ST1)
//          y = Vector.dense(FP64, A.nrows)
//          T.reduce_vector(out=y, accum=FP64.plus)
//          k = y.reduce_float()
//          return(1/k)*(3*A.plus_second(y) -2*T.plus_second(y) + y)

// TC3: in python:
//
//      def triangle_centrality3(A):
//          L = A.tril(-1)
//          T = A.plus_pair(A, mask=L, desc=descriptor.ST1)
//          T_T = T.T
//          y = T.reduce() + T_T.reduce()
//          k = y.reduce_float()
//          return (3 * A.plus_second(y) - (2 * (T.plus_second(y)
//                                           + T_T.plus_second(y))) + y) / k

// Note: the python TC3 forms T_T explicitly; but METHOD 3 below uses the
// descriptor to transpose T.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK           \
{                                   \
    GrB_free (&T) ;                 \
    GrB_free (&u) ;                 \
    GrB_free (&w) ;                 \
    GrB_free (&y) ;                 \
    GrB_free (&L) ;                 \
    GrB_free (&thunk) ;             \
}

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_FREE_WORK ;             \
    GrB_free (centrality) ;         \
}

#include "LG_internal.h"

// METHOD is 1, 15, 2, or 3 to select the above methods
#define METHOD 1

#if ( !LG_SUITESPARSE )
// methods 2 and 3 require SuiteSparse
#if ( METHOD >= 2 )
    #undef  METHOD
    #define METHOD 15
#endif
#endif

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
    GrB_Matrix T = NULL, L = NULL, A = NULL ;
    GrB_Vector y = NULL, u = NULL, w = NULL ;
    #if LG_SUITESPARSE
    GxB_Scalar thunk = NULL ;
    #else
    // not used, just to allow GrB_free to be done, in LAGraph_FREE_WORK
    GrB_Vector thunk = NULL ;
    #endif

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

    // FIXME: could remove any self-edges, if present; do this in the
    // non-expert version.
    // no self edges can be present
    LG_CHECK (G->ndiag != 0, -104, "G->ndiag must be zero") ;

// GxB_set (GxB_BURBLE, true) ;    // FIXME: remove

    //--------------------------------------------------------------------------
    // count triangles: T<A> = A*A' using the plus_pair semiring
    //--------------------------------------------------------------------------

    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_new (&T, GrB_FP64, n, n)) ;

#if (METHOD == 1 || METHOD == 15)

    //--------------------------------------------------------------------------
    // TC1, TC1.5: simplest method, requires that A has all entries equal to 1
    //--------------------------------------------------------------------------

//          T = A.mxm(A, mask=A)
//          y = T.reduce_vector()
//          k = y.reduce_float()
//          T = pattern of T        FIXME
//          return(1/k)*(3*(A @ y) - 2*(T @ y) + y)


    #if ( METHOD == 1)
        printf ("TC1: ")  ; // FIXME: remove
        // T<A> = A*A : method 1
        GrB_TRY (GrB_mxm (T, A, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A, A,
            NULL)) ;
    #else
        // this should be faster than METHOD 1
        printf ("TC1.5: ") ;   // FIXME: remove
        // T<A> = A*A' : method 1.5
        GrB_TRY (GrB_mxm (T, A, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A, A,
            GrB_DESC_T1)) ;
    #endif

    // y = sum (T), where y(i) = sum (T (i,:)) and y(i)=0 of T(i,:) is empty
    GrB_TRY (GrB_Vector_new (&y, GrB_FP64, n)) ;
    GrB_TRY (GrB_reduce (y, NULL, NULL, GrB_PLUS_MONOID_FP64, T, NULL)) ;

    // k = sum (y)
    double k = 0 ;
    GrB_TRY (GrB_reduce (&k, NULL, GrB_PLUS_MONOID_FP64, y, NULL)) ;

    // T = spones (T)
    GrB_TRY (GrB_assign (T, T, NULL, (double) 1, GrB_ALL, n, GrB_ALL, n,
        GrB_DESC_S)) ;

    // centrality = (3*A*y - 2*T*y + y) / k

    // w = T*y
    GrB_TRY (GrB_Vector_new (&w, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (w, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, T, y, NULL));

    // w = (-2)*w
    double minus_two = -2 ;
    GrB_TRY (GrB_apply (w, NULL, NULL, GrB_TIMES_FP64, minus_two, w, NULL)) ;

    // u = A*y
    GrB_TRY (GrB_Vector_new (&u, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (u, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A, y, NULL));

#elif (METHOD == 2)

    //--------------------------------------------------------------------------
    // TC2: using PLUS_PAIR semiring
    //--------------------------------------------------------------------------

    // only uses the pattern of A

    printf ("TC2: ")  ;     // FIXME: remove

    // T{A} = A*A'
    GrB_TRY (GrB_mxm (T, A, NULL, GxB_PLUS_PAIR_FP64, A, A, GrB_DESC_ST1)) ;

    // y = sum (T), where y(i) = sum (T (i,:)) and y(i)=0 of T(i,:) is empty
    GrB_TRY (GrB_Vector_new (&y, GrB_FP64, n)) ;
    GrB_TRY (GrB_assign (y, NULL, NULL, ((double) 0), GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_reduce (y, NULL, GrB_PLUS_FP64, GrB_PLUS_MONOID_FP64, T,
        NULL)) ;

    // k = sum (y)
    double k = 0 ;
    GrB_TRY (GrB_reduce (&k, NULL, GrB_PLUS_MONOID_FP64, y, NULL)) ;

    // centrality = (3*A*y - 2*T*y + y) / k

    // w = T*y
    GrB_TRY (GrB_Vector_new (&w, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (w, NULL, NULL, GxB_PLUS_SECOND_FP64, T, y, NULL)) ;

    // w = (-2)*w
    double minus_two = -2 ;
    GrB_TRY (GrB_apply (w, NULL, NULL, GrB_TIMES_FP64, minus_two, w, NULL)) ;

    // u = A*y
    GrB_TRY (GrB_Vector_new (&u, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (u, NULL, NULL, GxB_PLUS_SECOND_FP64, A, y, NULL)) ;

#elif (METHOD == 3)

    //--------------------------------------------------------------------------
    // TC3: using tril.  This is the fastest method.
    //--------------------------------------------------------------------------

    // only uses the pattern of A

    printf ("TC3: ")  ;     // FIXME: remove

    GrB_TRY (GrB_Matrix_new (&L, GrB_FP64, n, n)) ;

    // L = tril (A,-1)
    GrB_TRY (GxB_Scalar_new (&thunk, GrB_INT64)) ;
    GrB_TRY (GxB_Scalar_setElement (thunk, (int64_t) (-1))) ;
    GrB_TRY (GxB_select (L, NULL, NULL, GxB_TRIL, A, thunk, NULL)) ;
    GrB_TRY (GrB_free (&thunk)) ;

    // T{L}= A*A'
    GrB_TRY (GrB_mxm (T, L, NULL, GxB_PLUS_PAIR_FP64, A, A, GrB_DESC_ST1)) ;
    GrB_TRY (GrB_free (&L)) ;

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

    // w = (-2)*w
    double minus_two = -2 ;
    GrB_TRY (GrB_apply (w, NULL, NULL, GrB_TIMES_FP64, minus_two, w, NULL)) ;

    // u = A*y
    GrB_TRY (GrB_Vector_new (&u, GrB_FP64, n)) ;
    GrB_TRY (GrB_mxv (u, NULL, NULL, GxB_PLUS_SECOND_FP64, A, y, NULL)) ;

#endif

    //--------------------------------------------------------------------------
    // centrality = (3*u + w + y) / k for all 4 methods
    //--------------------------------------------------------------------------

    // centrality = 3*u
    GrB_TRY (GrB_Vector_new (centrality, GrB_FP64, n)) ;
    const double three = 3 ;
    GrB_TRY (GrB_apply (*centrality, NULL, NULL, GrB_TIMES_FP64, three, u,
        NULL)) ;

    // centrality += (w + y)
    GrB_TRY (GrB_eWiseAdd (*centrality, NULL, GrB_PLUS_FP64, GrB_PLUS_FP64,
        w, y, NULL)) ;

    // centrality = centrality / k
    GrB_TRY (GrB_apply (*centrality, NULL, NULL, GrB_TIMES_FP64,
        ((k == 0) ? 1.0 : (1.0/k)), *centrality, NULL)) ;

    // TODO: # of triangles is k/6, which could be returned as well
    printf ("# of triangles: %.32g\n", k/6) ;   // FIXME: remove

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;
// GxB_set (GxB_BURBLE, false) ;   // FIXME: remove
    return (0) ;
}

