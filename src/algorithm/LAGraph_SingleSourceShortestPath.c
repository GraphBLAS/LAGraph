//------------------------------------------------------------------------------
// LAGraph_SingleSourceShortestPath: single-source shortest path
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Single source shortest path with delta stepping.  Contributed by Jinhao
// Chen, Scott Kolodziej and Tim Davis, Texas A&M University.  Adapted from
// GraphBLAS Template Library (GBTL) by Scott McMillian and Tze Meng Low.

// U. Sridhar, M. Blanco, R. Mayuranath, D. G. Spampinato, T. M. Low, and
// S. McMillan, “Delta-Stepping SSSP: From Vertices and Edges to GraphBLAS
// Implementations,” in 2019 IEEE International Parallel and Distributed
// Processing Symposium Workshops (IPDPSW), 2019, pp. 241–250.
// https://ieeexplore.ieee.org/document/8778222/references
// https://arxiv.org/abs/1911.06895

// LAGraph_SingleSourceShortestPath computes the shortest path lengths from the
// specified source vertex to all other vertices in the graph.

// The parent vector is not computed; see LAGraph_BF_* instead.

// TODO: this method gets stuck in an infinite loop when there are negative-
// weight cycles in the graph.

#define LAGraph_FREE_WORK   \
{                           \
    GrB_free (&AL) ;        \
    GrB_free (&AH) ;        \
    GrB_free (&lBound) ;    \
    GrB_free (&uBound) ;    \
    GrB_free (&tmasked) ;   \
    GrB_free (&tReq) ;      \
    GrB_free (&tless) ;     \
    GrB_free (&s) ;         \
    GrB_free (&reach) ;     \
    GrB_free (&Empty) ;     \
}

#define LAGraph_FREE_ALL    \
{                           \
    LAGraph_FREE_WORK ;     \
    GrB_free (&t) ;         \
}

#include "LG_internal.h"

// TODO assert the input matrix has type GrB_INT32, or select different
// operators / semirings based on the matrix type.

int LAGraph_SingleSourceShortestPath    // returns 0 if successful, -1 if fail
(
    // output:
    GrB_Vector *path_length,    // path_length (i) is the length of the shortest
                                // path from the source vertex to vertex i
    // inputs:
    LAGraph_Graph G,
    GrB_Index source,           // source vertex
    int32_t delta,              // delta value for delta stepping
                                // TODO: use GrB_Scalar for delta
    // TODO: make this an enum, and add to LAGraph_Graph properties, and then
    // remove it from the inputs to this function
    //      case 0: A can have negative, zero, or positive entries
    //      case 1: A can have zero or positive entries
    //      case 2: A only has positive entries (see FIXME below)
    bool AIsAllPositive,       // A boolean indicating whether the entries of
                               // matrix A are all positive
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Scalar lBound = NULL ;  // the threshold for GrB_select
    GrB_Scalar uBound = NULL ;  // the threshold for GrB_select
    GrB_Matrix AL = NULL ;      // graph containing the light weight edges
    GrB_Matrix AH = NULL ;      // graph containing the heavy weight edges
    GrB_Vector t = NULL ;       // tentative shortest path length
    GrB_Vector tmasked = NULL ;
    GrB_Vector tReq = NULL ;
    GrB_Vector tless = NULL ;
    GrB_Vector s = NULL ;
    GrB_Vector reach = NULL ;
    GrB_Vector Empty = NULL ;

    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;
    LG_CHECK (path_length == NULL, -1, "path_length parameter is NULL") ;
    (*path_length) = NULL ;

    GrB_Matrix A = G->A ;
    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;

    LG_CHECK (source >= n || source < 0, -1, "invalid source node") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_Scalar_new (&lBound, GrB_INT32)) ;     // TODO: any type
    GrB_TRY (GrB_Scalar_new (&uBound, GrB_INT32)) ;     // TODO: any type
    GrB_TRY (GrB_Scalar_setElement (lBound, delta)) ;

    GrB_TRY (GrB_Vector_new (&t, GrB_INT32, n)) ;       // TODO: any type
    GrB_TRY (GrB_Vector_new (&tmasked, GrB_INT32, n)) ; // TODO: any type
    GrB_TRY (GrB_Vector_new (&tReq, GrB_INT32, n)) ;    // TODO: any type
    GrB_TRY (GrB_Vector_new (&Empty, GrB_BOOL, n)) ;

    GrB_TRY (GrB_Vector_new (&tless, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&s, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&reach, GrB_BOOL, n)) ;

#if LG_SUITESPARSE
    // optional hints for SuiteSparse:GraphBLAS
    GrB_TRY (GxB_set (t, GxB_SPARSITY_CONTROL, GxB_BITMAP)) ;
    GrB_TRY (GxB_set (tmasked, GxB_SPARSITY_CONTROL, GxB_SPARSE)) ;
    GrB_TRY (GxB_set (tReq, GxB_SPARSITY_CONTROL, GxB_SPARSE)) ;
    GrB_TRY (GxB_set (tless, GxB_SPARSITY_CONTROL, GxB_SPARSE)) ;
    GrB_TRY (GxB_set (s, GxB_SPARSITY_CONTROL, GxB_SPARSE)) ;
    GrB_TRY (GxB_set (reach, GxB_SPARSITY_CONTROL, GxB_BITMAP)) ;
#endif

    // t (:) = infinity
    GrB_TRY (GrB_assign (t, NULL, NULL,
        (int32_t) INT32_MAX,        // TODO: any type
        GrB_ALL, n, NULL)) ;
    // t (src) = 0
    GrB_TRY (GrB_Vector_setElement (t, 0, source)) ;

    // reach (src) = true
    GrB_TRY (GrB_Vector_setElement (reach, true, source)) ;

    // Instead of using tmasked >= i*delta = 0 to find out how many left to be
    // optimized, tmasked can be directly set the same as t since there is only
    // one entry that satisfies the condition.
    GrB_TRY (GrB_Vector_setElement (tmasked, 0, source)) ;
    LAGraph_TRY (LAGraph_Vector_wait (tmasked, msg)) ;

    // s (src) = true
    GrB_TRY (GrB_Vector_setElement (s, true, source)) ;

    // AL = A .* (A <= delta)
    GrB_TRY (GrB_Matrix_new (&AL, GrB_INT32, n, n)) ;   // TODO: any type
    GrB_TRY (GrB_select (AL, NULL, NULL, GrB_VALUELE_INT32, A, lBound, NULL)) ;
    LAGraph_TRY (LAGraph_Matrix_wait (AL, msg)) ;

    // AH = A .* (A > delta)
    GrB_TRY (GrB_Matrix_new (&AH, GrB_INT32, n, n)) ;   // TODO: any type
    GrB_TRY (GrB_select (AH, NULL, NULL, GrB_VALUEGT_INT32, A, lBound, NULL)) ;
    LAGraph_TRY (LAGraph_Matrix_wait (AH, msg)) ;

    //--------------------------------------------------------------------------
    // while (t >= i*delta) not empty
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; ; i++)
    {

        //----------------------------------------------------------------------
        // tmasked = all entries in t<reach> that are less than (i+1)*delta
        //----------------------------------------------------------------------

        // tmasked<reach> = t
        GrB_TRY (GrB_Vector_clear (tmasked)) ;
        GrB_TRY (GrB_assign (tmasked, reach, NULL, t, GrB_ALL, n, NULL)) ;

        // tmasked = select (tmasked < (i+1)*delta)
        GrB_TRY (GrB_Scalar_setElement (uBound, (i+1) * delta)) ;
        GrB_TRY (GrB_select (tmasked, NULL, NULL, GrB_VALUELT_INT32, tmasked,
            uBound, NULL)) ;
        GrB_Index tmasked_nvals ;
        GrB_TRY (GrB_Vector_nvals (&tmasked_nvals, tmasked)) ;

        //----------------------------------------------------------------------
        // continue while the current bucket (tmasked) is not empty
        //----------------------------------------------------------------------

        while (tmasked_nvals > 0)
        {
            // tReq = AL' (min.+) tmasked
            GrB_TRY (GrB_vxm (tReq, NULL, NULL,
                GrB_MIN_PLUS_SEMIRING_INT32,        // TODO: any type
                tmasked, AL, NULL)) ;

            // s<tmasked> = true
            GrB_TRY (GrB_assign (s, tmasked, NULL, (bool) true, GrB_ALL, n,
                GrB_DESC_S)) ;

            // if nnz (tReq) == 0, no need to continue the rest of this loop
            GrB_Index tReq_nvals ;
            GrB_TRY (GrB_Vector_nvals (&tReq_nvals, tReq)) ;
            if (tReq_nvals == 0) break ;

            // TODO currently assuming all edges weights are nonzero, so we're
            // using a structural mask.  Explicit zeros in A would require a
            // valued mask.
            // tless<tReq> = tReq .< t
            GrB_TRY (GrB_Vector_clear (tless)) ;
            GrB_TRY (GrB_eWiseAdd (tless, tReq, NULL,
                GrB_LT_INT32,        // TODO: any type
                tReq, t, GrB_DESC_S  /* assumes all entries in A are > 0 */)) ;

            // remove explicit zeros from tless so it can be used as a
            // structural mask
            GrB_Index tless_nvals ;
            GrB_TRY (GrB_select (tless, NULL, NULL, GrB_VALUENE_INT32, tless,
                (int32_t) 0, NULL)) ;
            GrB_TRY (GrB_Vector_nvals (&tless_nvals, tless)) ;
            if (tless_nvals == 0) break ;

            // update reachable node list/mask
            GrB_TRY (GrB_assign (reach, tless, NULL, (bool) true, GrB_ALL, n,
                GrB_DESC_S)) ;

            // tmasked<tless> = select (i*delta <= tReq < (i+1)*delta)
            // since all entries of the 5 GAP graphs are known to be
            // positive, and the entries of tmasked are at least i*delta,
            // tReq = tmasked min.+ AL must be >= i*delta.
            // Therefore, there is no need to perform GrB_select with
            // GrB_VALUEGE_INT32 to find tmasked >= i*delta from tReq 
            GrB_TRY (GrB_Vector_clear (tmasked)) ;
            GrB_TRY (GrB_select (tmasked, tless, NULL, GrB_VALUELT_INT32,
                tReq, uBound, GrB_DESC_S /* GrB_DESC_RS */)) ;

            // For general graph with some negative weights:
            if (!AIsAllPositive)
            {
                GrB_TRY (GrB_Scalar_setElement (lBound, i * delta)) ;
                // tmasked = select entries in tmasked that are >= lBound
                GrB_TRY (GrB_select (tmasked, NULL, NULL, GrB_VALUEGE_INT32, 
                    tmasked, lBound, NULL)) ;
            }

            // t<tless> = tReq
            // TODO: use assign, not apply
            GrB_TRY (GrB_apply (t, tless, NULL,
                GrB_IDENTITY_INT32,     // TODO: any type
                tReq, GrB_DESC_S)) ;
            GrB_TRY (GrB_Vector_nvals (&tmasked_nvals, tmasked)) ;
        }

        // tmasked<s> = t
        GrB_TRY (GrB_assign (tmasked, s, NULL, t, GrB_ALL, n, GrB_DESC_RS)) ;

        // tReq = AH'*tmasked
        GrB_TRY (GrB_vxm (tReq, NULL, NULL,
            GrB_MIN_PLUS_SEMIRING_INT32,    // TODO: any type
            tmasked, AH, NULL)) ;

        // t = min (t, tReq)
        // When t is dense, it is best to get tless<tReq> = tReq .< t,
        // and use tless as mask to update t.
        GrB_TRY (GrB_Vector_clear (tless)) ;
        GrB_TRY (GrB_eWiseAdd (tless, tReq, NULL,
            GrB_LT_INT32,           // TODO: any type
            tReq, t, GrB_DESC_S)) ;
        // t<tless> = tReq
        // TODO: use assign, not apply
        GrB_TRY (GrB_apply (t, tless, NULL,
            GrB_IDENTITY_INT32,     // TODO: any type
            tReq, NULL)) ;

        //----------------------------------------------------------------------
        // find out how many left to be computed
        //----------------------------------------------------------------------

        // update reachable node list/mask
        GrB_TRY (GrB_assign (reach, tless, NULL, (bool) true, GrB_ALL, n,
            NULL)) ;

        // remove previous buckets
        // reach<s,struct> = Empty
        GrB_TRY (GrB_assign (reach, s, NULL, Empty, GrB_ALL, n, GrB_DESC_S)) ;
        GrB_Index nreach ;
        GrB_TRY (GrB_Vector_nvals (&nreach, reach)) ;
        if (nreach == 0) break ;

        GrB_TRY (GrB_Vector_clear (s)) ; // clear s for the next iteration
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*path_length) = t ;
    LAGraph_FREE_WORK ;
    return (0) ;
}

