//------------------------------------------------------------------------------
// LAGraph_sssp12: single-source shortest path
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2020 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------
// LAGraph_sssp13: Single source shortest path with delta stepping.
// Contributed by Jinhao Chen, Scott Kolodziej and Tim Davis, Texas A&M
// University.
// Adapted from GraphBLAS Template Library (GBTL) by Scott McMillian and Tze
// Meng Low.
//------------------------------------------------------------------------------
// U. Sridhar, M. Blanco, R. Mayuranath, D. G. Spampinato, T. M. Low, and
// S. McMillan, “Delta-Stepping SSSP: From Vertices and Edges to GraphBLAS
// Implementations,” in 2019 IEEE International Parallel and Distributed
// Processing Symposium Workshops (IPDPSW), 2019, pp. 241–250.
// https://ieeexplore.ieee.org/document/8778222/references
// https://arxiv.org/abs/1911.06895
//------------------------------------------------------------------------------
// LAGraph_sssp computes the shortest path lengths from the specified source
// vertex to all other vertices in the graph.
//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define TIMING 0

#define LAGRAPH_FREE_ALL                    \
{                                           \
    GrB_free(&AL);                          \
    GrB_free(&AH);                          \
    GrB_free(&ALT);                         \
    GrB_free(&AHT);                         \
    GrB_free(&lBound);                      \
    GrB_free(&uBound);                      \
    GrB_free(&Inf);                         \
    GrB_free(&t);                           \
    GrB_free(&tmasked);                     \
    GrB_free(&tReq);                        \
    GrB_free(&tless);                       \
    GrB_free(&s);                           \
    GrB_free(&reach);                       \
}

// TODO assert the input matrix has type GrB_INT32
// TODO try different delta (for the GAP matrices)

GrB_Info LAGraph_sssp13         // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in
                               // semiring (INT32)
    GrB_Matrix AT,             // AT = A'
    GrB_Index source,          // source vertex from which to compute
                               // shortest paths
    int32_t delta,             // delta value for delta stepping
    // TODO: make this an enum:
    //      case 0: A can have negative, zero, or positive entries
    //      case 1: A can have zero or positive entries
    //      case 2: A only has positive entries (see FIXME below)
    bool AIsAllPositive        // A boolean indicating whether the entries of
                               // matrix A are all positive
)
{

    GrB_Info info;

    (*path_length) = NULL;
    GrB_Index nrows, ncols, n = 0; // graph info
    GrB_Index ignore, tmasked_nvals = 0, s_nvals;

    GxB_Scalar lBound = NULL; // the threshold for GxB_select
    GxB_Scalar uBound = NULL; // the threshold for GxB_select
    GxB_Scalar Inf = NULL; // the max value allowed in t, i.e., INT32_MAX

    GrB_Matrix AL = NULL;     // graph containing the light weight edges
    GrB_Matrix AH = NULL;     // graph containing the heavy weight edges

    GrB_Matrix ALT = NULL;    // ALT = AL'
    GrB_Matrix AHT = NULL;    // AHT = AH'

    // GrB_INT32 vectors
    GrB_Vector t = NULL;    // tentative shortest path length
    GrB_Vector tmasked = NULL;
    GrB_Vector tReq = NULL;

    // GrB_BOOL vectors
    GrB_Vector tless = NULL;
    GrB_Vector s = NULL;
    GrB_Vector reach = NULL;

    // decision on the first operation among t>=i*delta or t<Inf
    // initially, there tends to have more Inf, so we do t<Inf first to
    // get sparser result
    bool do_LT_first = true;

    if (A == NULL || path_length == NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    // Get dimensions
    LAGr_Matrix_nrows (&nrows, A) ;
    LAGr_Matrix_ncols (&ncols, A) ;

    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE) ;
    }

    n = nrows;

    if (source >= n || source < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex", GrB_INVALID_VALUE) ;
    }

    LAGr_Scalar_new(&lBound, GrB_INT32) ;
    LAGr_Scalar_new(&uBound, GrB_INT32) ;
    LAGr_Scalar_new(&Inf, GrB_INT32) ;
    LAGr_Scalar_setElement (lBound, delta) ;
    LAGr_Scalar_setElement (Inf, INT32_MAX) ;

    // Create the workspace vectors
    LAGr_Vector_new(&t, GrB_INT32, n);
    LAGr_Vector_new(&tmasked, GrB_INT32, n);
    LAGr_Vector_new(&tReq, GrB_INT32, n);

    LAGr_Vector_new(&tless, GrB_BOOL, n);
    LAGr_Vector_new(&s, GrB_BOOL, n);
    LAGr_Vector_new(&reach, GrB_BOOL, n);

    // t = infinity, t[src] = 0
    LAGr_assign(t, NULL, NULL, INT32_MAX, GrB_ALL, n, NULL);
    LAGr_Vector_setElement(t, 0, source);

    // reach = false, reach[src]= true
    LAGr_assign(reach, NULL, NULL, false, GrB_ALL, n, NULL);
    LAGr_Vector_setElement(reach, true, source);

    // if there is any node that can be further relaxed
    bool remain = true;

    // TODO: computing AL and AH could be done with one pass,
    // but it would require import/export.

    // AL = A .* (A <= delta) with lBound = delta
    LAGr_Matrix_new(&AL, GrB_INT32, n, n);
    LAGr_select(AL, NULL, NULL, GxB_LE_THUNK, A, lBound, NULL) ;

    // AL = AT'
    LAGr_Matrix_new(&ALT, GrB_INT32, n, n);
    LAGr_select(ALT, NULL, NULL, GxB_LE_THUNK, AT, lBound, NULL) ;

    // AH = A .* (A > delta) with lBound = delta
    LAGr_Matrix_new(&AH, GrB_INT32, n, n);
    LAGr_select(AH, NULL, NULL, GxB_GT_THUNK, A, lBound, NULL) ;

    // AHT = AH'
    LAGr_Matrix_new(&AHT, GrB_INT32, n, n);
    LAGr_select(AHT, NULL, NULL, GxB_GT_THUNK, AT, lBound, NULL) ;

    GrB_Index nvals = 0;
    LAGr_Matrix_nvals(&nvals, A);
    LAGr_Matrix_nvals(&nvals, AL);
    LAGr_Matrix_nvals(&nvals, AH);

    int32_t i = 0;

    // Instead of using tmasked >= i*delta = 0 to find out how many left to be
    // optimized, tmasked can be directly set the same as t since there is only
    // one entry that satisfies the condition
    // Furthermore, set s[src] = true to get the correct result from
    // GxB_PAIR_BOOL in the first loop.
    //LAGr_Scalar_setElement (lBound, i * delta);
    LAGr_Vector_setElement(tmasked, 0, source);
    LAGr_Vector_setElement(s, true, source);

    LAGr_Vector_nvals(&tmasked_nvals, tmasked);

    //--------------------------------------------------------------------------
    // while (t >= i*delta) not empty
    //--------------------------------------------------------------------------

    while (remain)
    {
        // tmasked = select (tmasked < (i+1)*delta)
        LAGr_Vector_clear (tmasked) ;
        LAGr_Scalar_setElement (uBound, (i+1) * delta);
        LAGr_assign (tmasked, reach, NULL, t, GrB_ALL, n, NULL) ;
        LAGr_select(tmasked, NULL, NULL, GxB_LT_THUNK, tmasked, uBound, NULL);

        LAGr_Vector_nvals(&tmasked_nvals, tmasked);

        //------------------------------------------------------------------
        // continue while the current bucket B[i] is not empty
        //------------------------------------------------------------------

        while (tmasked_nvals > 0)
        {
            // tReq = AL' (min.+) tmasked
            // int sparsity ;
            //if (sparsity == GxB_BITMAP || sparsity == GxB_FULL)
            if (tmasked_nvals > n/10)
            {
                // pull
                GxB_set (tmasked, GxB_SPARSITY_CONTROL, GxB_BITMAP) ;
                LAGr_mxv (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, ALT, tmasked, NULL) ;
            }
            else
            {
                // push
                GxB_set (tmasked, GxB_SPARSITY_CONTROL, GxB_SPARSE) ;
                LAGr_vxm (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, tmasked, AL, NULL) ;
            }

            // Even though GrB_assign is faster than eWiseAdd here, the
            // total time of getting s with assign and tmasked = (s.*t)
            // later is longer than using eWiseAdd here
            // Time taken here and time to get tmasked =(s.*t) is commented
            // below for kron matrix
            // s = (s | pattern of tmasked)
            LAGr_eWiseAdd (s, NULL, NULL, GxB_PAIR_BOOL, s, tmasked, NULL);

            // if nnz(tReq) == 0, no need to continue the rest of this loop
            GrB_Index tReq_nvals ;
            LAGr_Vector_nvals(&tReq_nvals, tReq);
            if (tReq_nvals == 0) { break; }

            // tless<tReq> = tReq .< t

            // TODO currently assuming all edges weights > 0 (FIXME), so 
            // we're using a structural mask.  Explicit zeros in A would
            // require a valued mask.

            LAGr_Vector_clear (tless) ;
            LAGr_eWiseAdd (tless, tReq, NULL, GrB_LT_INT32, tReq, t,
                GrB_DESC_S  /* assumes all entries in A are > 0 */) ;

            // remove explicit zeros from tless so it can be used as a
            // structural mask
            GrB_Index tless_nvals ;
            LAGr_select (tless, NULL, NULL, GxB_NONZERO, tless, NULL, NULL);
            LAGr_Vector_nvals (&tless_nvals, tless) ;
            if (tless_nvals == 0) { break ; }

            // update reachable node list/mask
            LAGr_assign (reach, tless, NULL, true, GrB_ALL, n, GrB_DESC_S) ;

            // tmasked<tless> = select (i*delta <= tReq < (i+1)*delta)
            // since all entries of the 5 GAP graphs are known to be
            // positive, and the entries of tmasked are at least i*delta,
            // tReq = tmasked min.+ AL must be >= i*delta.
            // Therefore, there is no need to perform GxB_select with
            // GxB_GE_THUNK to find tmasked >= i*delta from tReq 
            LAGr_Vector_clear (tmasked) ;
            LAGr_select (tmasked, tless, NULL, GxB_LT_THUNK,
                tReq, uBound, GrB_DESC_S /* GrB_DESC_RS */) ;
            // For general graph with negative weight, the following needs
            // to be done
            if (!AIsAllPositive)
            {
                LAGr_Scalar_setElement (lBound, i * delta);
                // tmasked = select entries in tmasked that are >= lBound
                LAGr_select (tmasked, NULL, NULL, GxB_GE_THUNK, 
                    tmasked, lBound, NULL) ;
            }

            // t<tless> = tReq
            LAGr_apply (t, tless, NULL, GrB_IDENTITY_INT32, tReq,
                GrB_DESC_S) ;

            LAGr_Vector_nvals(&tmasked_nvals, tmasked);
        }

        // tmasked<s> = t
        LAGr_assign (tmasked, s, NULL, t, GrB_ALL, n, GrB_DESC_RS) ;

        // tReq = AH'*tmasked
        // LAGr_vxm (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, tmasked, AH, NULL) ;
        // int sparsity ;
        // GxB_get (tmasked, GxB_SPARSITY_STATUS, &sparsity) ;
        // if (sparsity == GxB_BITMAP || sparsity == GxB_FULL)
        LAGr_Vector_nvals(&tmasked_nvals, tmasked);
        if (tmasked_nvals > n/10)
        {
            // pull
            LAGr_mxv (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, AHT, tmasked, NULL) ;
        }
        else
        {
            // push
            LAGr_vxm (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, tmasked, AH, NULL) ;
        }

        // t = min(t, tReq)
        // When t is dense, it is best to get tless<tReq> = tReq .< t,
        // and use tless as mask to update t.
        LAGr_Vector_clear (tless) ;
        LAGr_eWiseAdd (tless, tReq, NULL, GrB_LT_INT32, tReq, t,
            GrB_DESC_S);
        LAGr_apply(t, tless, NULL, GrB_IDENTITY_INT32, tReq, NULL);

        //------------------------------------------------------------------
        // find out how many left to be computed
        //------------------------------------------------------------------

        // update reachable node list/mask
        LAGr_assign (reach, tless, NULL, true, GrB_ALL, n, NULL) ;

        // remove previous buckets
        LAGr_assign (reach, s, NULL, false, GrB_ALL, n, GrB_DESC_S) ;
        LAGr_reduce (&remain, NULL, GxB_LOR_BOOL_MONOID, reach, NULL) ;

        LAGr_Vector_clear (s) ; // clear s for the next loop
        ++i;
    }

    // result = t
    *path_length = t;
    t = NULL;
    LAGRAPH_FREE_ALL ;
    return GrB_SUCCESS;
}

