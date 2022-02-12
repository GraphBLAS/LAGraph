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
// S. McMillan, "Delta-Stepping SSSP: From Vertices and Edges to GraphBLAS
// Implementations," in 2019 IEEE International Parallel and Distributed
// Processing Symposium Workshops (IPDPSW), 2019, pp. 241–250.
// https://ieeexplore.ieee.org/document/8778222/references
// https://arxiv.org/abs/1911.06895

// LAGraph_SingleSourceShortestPath computes the shortest path lengths from the
// specified source vertex to all other vertices in the graph.

// The parent vector is not computed; see LAGraph_BF_* instead.

// NOTE: this method gets stuck in an infinite loop when there are negative-
// weight cycles in the graph.

// TODO: add AIsAllPositive or related as a G->property.  But then this would
// not be a Basic method.  The Advanced method would require the AIsAllPostive
// cached property, and it would require Delta to be provided.

// TODO: Should a Basic method pick delta automatically?  The Basic method
// would compute the cached property AIsAllPositive (or related), and then
// it would also try to set Delta.  What should Delta be for an arbitrary
// graph, of type int32, int64, uint32, uint64, float, or double?  Perhaps
// equal some multiple (30?) of the max edge weight?  Unsure.

#define LG_FREE_WORK        \
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

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&t) ;         \
}

#include "LG_internal.h"

#define setelement(s, k)                                                      \
{                                                                             \
    switch (tcode)                                                            \
    {                                                                         \
        default:                                                              \
        case 0 : GrB_Scalar_setElement_INT32  (s, k * delta_int32 ) ; break ; \
        case 1 : GrB_Scalar_setElement_INT64  (s, k * delta_int64 ) ; break ; \
        case 2 : GrB_Scalar_setElement_UINT32 (s, k * delta_uint32) ; break ; \
        case 3 : GrB_Scalar_setElement_UINT64 (s, k * delta_uint64) ; break ; \
        case 4 : GrB_Scalar_setElement_FP32   (s, k * delta_fp32  ) ; break ; \
        case 5 : GrB_Scalar_setElement_FP64   (s, k * delta_fp64  ) ; break ; \
    }                                                                         \
    /* printf ("setelement: k %ld\n", k) */ ; \
    /* GxB_print (s, 3) */ ; \
}

int LAGraph_SingleSourceShortestPath
(
    // output:
    GrB_Vector *path_length,    // path_length (i) is the length of the shortest
                                // path from the source vertex to vertex i
    // input:
    LAGraph_Graph G,
    GrB_Index source,           // source vertex
    GrB_Scalar Delta,           // delta value for delta stepping
    // TODO: make this an enum, and add to LAGraph_Graph properties, and then
    // remove it from the inputs to this function
    //      case 0: A can have negative, zero, or positive entries
    //      case 1: A can have zero or positive entries
    //      case 2: A only has positive entries
    // TODO: add AIsAllPositive to G->A_is_something...
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

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (path_length != NULL && Delta != NULL, GrB_NULL_POINTER) ;
    GrB_Index nvals ;
    LG_TRY (GrB_Scalar_nvals (&nvals, Delta)) ;
    LG_ASSERT_MSG (nvals == 1, GrB_EMPTY_OBJECT, "delta is missing") ;
    (*path_length) = NULL ;

    GrB_Matrix A = G->A ;
    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    LG_ASSERT_MSG (source < n, GrB_INVALID_INDEX, "invalid source node") ;

    // GxB_print (Delta, 3) ;
    // GxB_print (A, 3) ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // get the type of the A matrix
    GrB_Type etype ;
    char typename [LAGRAPH_MAX_NAME_LEN] ;
    LG_TRY (LAGraph_Matrix_TypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&etype, typename, msg)) ;

    GrB_TRY (GrB_Scalar_new (&lBound, etype)) ;
    GrB_TRY (GrB_Scalar_new (&uBound, etype)) ;
    GrB_TRY (GrB_Vector_new (&t, etype, n)) ;
    GrB_TRY (GrB_Vector_new (&tmasked, etype, n)) ;
    GrB_TRY (GrB_Vector_new (&tReq, etype, n)) ;
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

    // select the operators, and set t (:) = infinity
    GrB_IndexUnaryOp ne, le, ge, lt, gt ;
    GrB_BinaryOp less_than ;
    GrB_Semiring min_plus ;
    int tcode ;
    int32_t  delta_int32  ;
    int64_t  delta_int64  ;
    uint32_t delta_uint32 ;
    uint64_t delta_uint64 ;
    float    delta_fp32   ;
    double   delta_fp64   ;

    // GxB_print (etype, 3) ;

    if (etype == GrB_INT32)
    {
        GrB_TRY (GrB_Scalar_extractElement (&delta_int32, Delta)) ;
        GrB_TRY (GrB_assign (t, NULL, NULL, (int32_t) INT32_MAX,
            GrB_ALL, n, NULL)) ;
        ne = GrB_VALUENE_INT32 ;
        le = GrB_VALUELE_INT32 ;
        ge = GrB_VALUEGE_INT32 ;
        lt = GrB_VALUELT_INT32 ;
        gt = GrB_VALUEGT_INT32 ;
        less_than = GrB_LT_INT32 ;
        min_plus = GrB_MIN_PLUS_SEMIRING_INT32 ;
        tcode = 0 ;
    }
    else if (etype == GrB_INT64)
    {
        GrB_TRY (GrB_Scalar_extractElement (&delta_int64, Delta)) ;
        GrB_TRY (GrB_assign (t, NULL, NULL, (int64_t) INT64_MAX,
            GrB_ALL, n, NULL)) ;
        ne = GrB_VALUENE_INT64 ;
        le = GrB_VALUELE_INT64 ;
        ge = GrB_VALUEGE_INT64 ;
        lt = GrB_VALUELT_INT64 ;
        gt = GrB_VALUEGT_INT64 ;
        less_than = GrB_LT_INT64 ;
        min_plus = GrB_MIN_PLUS_SEMIRING_INT64 ;
        tcode = 1 ;
    }
    else if (etype == GrB_UINT32)
    {
        GrB_TRY (GrB_Scalar_extractElement (&delta_uint32, Delta)) ;
        GrB_TRY (GrB_assign (t, NULL, NULL, (uint32_t) UINT32_MAX,
            GrB_ALL, n, NULL)) ;
        ne = GrB_VALUENE_UINT32 ;
        le = GrB_VALUELE_UINT32 ;
        ge = GrB_VALUEGE_UINT32 ;
        lt = GrB_VALUELT_UINT32 ;
        gt = GrB_VALUEGT_UINT32 ;
        less_than = GrB_LT_UINT32 ;
        min_plus = GrB_MIN_PLUS_SEMIRING_UINT32 ;
        AIsAllPositive = true ;
        tcode = 2 ;
    }
    else if (etype == GrB_UINT64)
    {
        GrB_TRY (GrB_Scalar_extractElement (&delta_uint64, Delta)) ;
        GrB_TRY (GrB_assign (t, NULL, NULL, (uint64_t) UINT64_MAX,
            GrB_ALL, n, NULL)) ;
        ne = GrB_VALUENE_UINT64 ;
        le = GrB_VALUELE_UINT64 ;
        ge = GrB_VALUEGE_UINT64 ;
        lt = GrB_VALUELT_UINT64 ;
        gt = GrB_VALUEGT_UINT64 ;
        less_than = GrB_LT_UINT64 ;
        min_plus = GrB_MIN_PLUS_SEMIRING_UINT64 ;
        AIsAllPositive = true ;
        tcode = 3 ;
    }
    else if (etype == GrB_FP32)
    {
        GrB_TRY (GrB_Scalar_extractElement (&delta_fp32, Delta)) ;
        GrB_TRY (GrB_assign (t, NULL, NULL, (float) INFINITY,
            GrB_ALL, n, NULL)) ;
        ne = GrB_VALUENE_FP32 ;
        le = GrB_VALUELE_FP32 ;
        ge = GrB_VALUEGE_FP32 ;
        lt = GrB_VALUELT_FP32 ;
        gt = GrB_VALUEGT_FP32 ;
        less_than = GrB_LT_FP32 ;
        min_plus = GrB_MIN_PLUS_SEMIRING_FP32 ;
        tcode = 4 ;
    }
    else if (etype == GrB_FP64)
    {
        GrB_TRY (GrB_Scalar_extractElement (&delta_fp64, Delta)) ;
        GrB_TRY (GrB_assign (t, NULL, NULL, (double) INFINITY,
            GrB_ALL, n, NULL)) ;
        ne = GrB_VALUENE_FP64 ;
        le = GrB_VALUELE_FP64 ;
        ge = GrB_VALUEGE_FP64 ;
        lt = GrB_VALUELT_FP64 ;
        gt = GrB_VALUEGT_FP64 ;
        less_than = GrB_LT_FP64 ;
        min_plus = GrB_MIN_PLUS_SEMIRING_FP64 ;
        tcode = 5 ;
    }
    else
    {
        LG_ASSERT_MSG (false, GrB_NOT_IMPLEMENTED, "type not supported") ;
    }

    // t (src) = 0
    GrB_TRY (GrB_Vector_setElement (t, 0, source)) ;

    // reach (src) = true
    GrB_TRY (GrB_Vector_setElement (reach, true, source)) ;

    // Instead of using tmasked >= i*delta = 0 to find out how many left to be
    // optimized, tmasked can be directly set the same as t since there is only
    // one entry that satisfies the condition.
    GrB_TRY (GrB_Vector_setElement (tmasked, 0, source)) ;
    GrB_TRY (GrB_wait (tmasked, GrB_MATERIALIZE)) ;

    // s (src) = true
    GrB_TRY (GrB_Vector_setElement (s, true, source)) ;

    // AL = A .* (A <= delta)
    GrB_TRY (GrB_Matrix_new (&AL, etype, n, n)) ;
    GrB_TRY (GrB_select (AL, NULL, NULL, le, A, Delta, NULL)) ;
    GrB_TRY (GrB_wait (AL, GrB_MATERIALIZE)) ;

    // AH = A .* (A > delta)
    GrB_TRY (GrB_Matrix_new (&AH, etype, n, n)) ;
    GrB_TRY (GrB_select (AH, NULL, NULL, gt, A, Delta, NULL)) ;
    GrB_TRY (GrB_wait (AH, GrB_MATERIALIZE)) ;

    //--------------------------------------------------------------------------
    // while (t >= i*delta) not empty
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; ; i++)
    {

        //----------------------------------------------------------------------
        // tmasked = all entries in t<reach> that are less than (i+1)*delta
        //----------------------------------------------------------------------

        // tmasked<reach> = t
        GrB_TRY (GrB_Vector_clear (tmasked)) ;  // TODO: delete this
        GrB_TRY (GrB_assign (tmasked, reach, NULL, t, GrB_ALL, n,
            NULL)) ;    // TODO: use GrB_DESC_RS

        // tmasked = select (tmasked < (i+1)*delta)
        setelement (uBound, (i+1)) ;        // uBound = (i+1) * Delta
        GrB_TRY (GrB_select (tmasked, NULL, NULL, lt, tmasked, uBound, NULL)) ;
        GrB_Index tmasked_nvals ;
        GrB_TRY (GrB_Vector_nvals (&tmasked_nvals, tmasked)) ;

        //----------------------------------------------------------------------
        // continue while the current bucket (tmasked) is not empty
        //----------------------------------------------------------------------

        while (tmasked_nvals > 0)
        {
            // tReq = AL'*tmasked using the min_plus semiring
            GrB_TRY (GrB_vxm (tReq, NULL, NULL, min_plus, tmasked, AL, NULL)) ;

            // s<struct(tmasked)> = true
            GrB_TRY (GrB_assign (s, tmasked, NULL, (bool) true, GrB_ALL, n,
                GrB_DESC_S)) ;

            // if nnz (tReq) == 0, no need to continue the rest of this loop
            GrB_Index tReq_nvals ;
            GrB_TRY (GrB_Vector_nvals (&tReq_nvals, tReq)) ;
            if (tReq_nvals == 0) break ;

            // tless<tReq> = tReq .< t
            // TODO: try eWiseMult with no mask
            GrB_TRY (GrB_Vector_clear (tless)) ;    // TODO: delete this
            GrB_TRY (GrB_eWiseAdd (tless, tReq, NULL, less_than, tReq, t,
                GrB_DESC_S)) ;  // TODO: use GrB_DESC_RS

            // remove explicit zeros from tless so it can be used as a
            // structural mask
            GrB_Index tless_nvals ;
            GrB_TRY (GrB_select (tless, NULL, NULL, ne, tless, (int32_t) 0,
                NULL)) ;
            GrB_TRY (GrB_Vector_nvals (&tless_nvals, tless)) ;
            if (tless_nvals == 0) break ;

            // update reachable node list/mask
            // reach<struct(tless)> = true
            GrB_TRY (GrB_assign (reach, tless, NULL, (bool) true, GrB_ALL, n,
                GrB_DESC_S)) ;

            // tmasked<tless> = select (i*delta <= tReq < (i+1)*delta).
            // If all entries of the graph is known to be positive, and the
            // entries of tmasked are at least i*delta, tReq = tmasked min.+ AL
            // must be >= i*delta.  Therefore, there is no need to perform
            // GrB_select with ge to find tmasked >= i*delta from tReq
            GrB_TRY (GrB_Vector_clear (tmasked)) ;  // TODO: delete this
            GrB_TRY (GrB_select (tmasked, tless, NULL, lt, tReq, uBound,
                GrB_DESC_S)) ;  // TODO: use GrB_DESC_RS

            // For general graph with some negative weights:
            if (!AIsAllPositive)
            {
                setelement (lBound, (i)) ;  // lBound = i * Delta
                // tmasked = select entries in tmasked that are >= lBound
                GrB_TRY (GrB_select (tmasked, NULL, NULL, ge, tmasked, lBound,
                    NULL)) ;
            }

            // t<struct(tless)> = tReq
            GrB_TRY (GrB_assign (t, tless, NULL, tReq, GrB_ALL, n, GrB_DESC_S));
            GrB_TRY (GrB_Vector_nvals (&tmasked_nvals, tmasked)) ;
        }

        // tmasked<s> = t
        GrB_TRY (GrB_assign (tmasked, s, NULL, t, GrB_ALL, n, GrB_DESC_RS)) ;

        // tReq = AH'*tmasked using the min_plus semiring
        GrB_TRY (GrB_vxm (tReq, NULL, NULL, min_plus, tmasked, AH, NULL)) ;

        // t = min (t, tReq)
        // When t is dense, it is best to get tless<tReq> = tReq .< t,
        // and use tless as mask to update t.
        GrB_TRY (GrB_Vector_clear (tless)) ;    // TODO: delete this
        GrB_TRY (GrB_eWiseAdd (tless, tReq, NULL, less_than, tReq, t,
            GrB_DESC_S)) ;  // GrB_DESC_RS, or try eWiseMult with no mask
        // t<tless> = tReq
        GrB_TRY (GrB_assign (t, tless, NULL, tReq, GrB_ALL, n, NULL)) ;

        //----------------------------------------------------------------------
        // find out how many left to be computed
        //----------------------------------------------------------------------

        // update reachable node list/mask
        // reach<tless> = true
        GrB_TRY (GrB_assign (reach, tless, NULL, (bool) true, GrB_ALL, n,
            NULL)) ;

        // remove previous buckets
        // reach<struct(s)> = Empty
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
    LG_FREE_WORK ;
    // GxB_print (t, 3) ;
    return (GrB_SUCCESS) ;
}

