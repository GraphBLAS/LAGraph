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
// LAGraph_sssp12: Single source shortest path with delta stepping.
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

GrB_Info LAGraph_sssp12c        // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in
                               // semiring (INT32)
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
    // int32_t print_lvl = 0; // change to 2 to show calculation step results

    (*path_length) = NULL;
    GrB_Index nrows, ncols, n = 0; // graph info
    GrB_Index ignore, tmasked_nvals = 0, s_nvals;

    GxB_Scalar lBound = NULL; // the threshold for GxB_select
    GxB_Scalar uBound = NULL; // the threshold for GxB_select
    GxB_Scalar Inf = NULL; // the max value allowed in t, i.e., INT32_MAX

    GrB_Matrix AL = NULL;     // graph containing the light weight edges
    GrB_Matrix AH = NULL;     // graph containing the heavy weight edges

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

    // double total_time1 = 0;
    // double total_time2 = 0;
    // double total_time3 = 0;
    // double total_time4 = 0;
    // double total_time5 = 0;
    // double total_time6 = 0;
    // double total_time7 = 0;
    // double total_time9 = 0;
    // double total_time10= 0;
    // double tic1[2], tic[2];
    // LAGraph_tic(tic1);
    // LAGraph_tic(tic);

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
    // GxB_print(AL, print_lvl);

    // AH = A .* (A > delta) with lBound = delta
    LAGr_Matrix_new(&AH, GrB_INT32, n, n);
    LAGr_select(AH, NULL, NULL, GxB_GT_THUNK, A, lBound, NULL) ;

    GrB_Index nvals = 0;
    LAGr_Matrix_nvals(&nvals, A);
    // printf("A has nvals = %"PRIu64"\n", nvals);
    LAGr_Matrix_nvals(&nvals, AL);
    // printf("AL has nvals = %"PRIu64"\n", nvals);
    LAGr_Matrix_nvals(&nvals, AH);
    // printf("AH has nvals = %"PRIu64"\n", nvals);
    // GxB_print(AH, print_lvl);

    int32_t i = 0;

    // Instead of using tmasked >= i*delta = 0 to find out how many left to be
    // optimized, tmasked can be directly set the same as t since there is only
    // one entry that satisfies the condition
    // Furthermore, set s[src] = true to get the correct result from
    // GxB_PAIR_BOOL in the first loop.
    //LAGr_Scalar_setElement (lBound, i * delta);
    LAGr_Vector_setElement(tmasked, 0, source);
    LAGr_Vector_setElement(s, true, source);
    // GxB_print(tmasked, print_lvl);

    LAGr_Vector_nvals(&tmasked_nvals, tmasked);

//    if (print_lvl > 0)
//    {
//        fprintf (stderr, "outter tmasked has %ld nnz\n",tmasked_nvals);
//    }
    // double t_pre = LAGraph_toc(tic);
    //printf("pre-handling time %12.6g sec\n", t1);
    // double t1 ;

    //--------------------------------------------------------------------------
    // while (t >= i*delta) not empty
    //--------------------------------------------------------------------------

    while (remain)
    {
        // printf ("\n============================= outer: %d\n", i) ;
        // tmasked = select (tmasked < (i+1)*delta)
        // LAGraph_tic (tic);
        LAGr_Vector_clear (tmasked) ;
        LAGr_Scalar_setElement (uBound, (i+1) * delta);
        LAGr_assign (tmasked, reach, NULL, t, GrB_ALL, n, NULL) ;
        LAGr_select(tmasked, NULL, NULL, GxB_LT_THUNK, tmasked, uBound, NULL);
        // t1 = LAGraph_toc(tic);
        // total_time1 += t1;

        LAGr_Vector_nvals(&tmasked_nvals, tmasked);

        //------------------------------------------------------------------
        // continue while the current bucket B[i] is not empty
        //------------------------------------------------------------------

        while (tmasked_nvals > 0)
        {
            //printf ("\n=============== inner tmasked has %ld nnz\n",tmasked_nvals);
            // tReq = AL' (min.+) tmasked
            // LAGraph_tic (tic);
            LAGr_vxm (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, tmasked, AL, NULL) ;
            // t1 = LAGraph_toc(tic);
            // total_time2 += t1;
            //GxB_print(tReq, 2);

            // Even though GrB_assign is faster than eWiseAdd here, the
            // total time of getting s with assign and tmasked = (s.*t)
            // later is longer than using eWiseAdd here
            // Time taken here and time to get tmasked =(s.*t) is commented
            // below for kron matrix
            // s = (s | pattern of tmasked)
            // LAGraph_tic (tic);

            //printf("-------------------------------------------------\n");
            //GxB_print(s, 2);
            //GxB_print(tmasked, 2);
            // best:0.35+1.17=1.52sec
            LAGr_eWiseAdd (s, NULL, NULL, GxB_PAIR_BOOL, s, tmasked, NULL);
            // 2nd best:0.50+1.22=1.72sec
            //LAGr_apply(s, NULL, GrB_LOR, GxB_ONE_BOOL, tmasked, NULL);
            // worse:0.14+2.27=2.41sec
            //LAGr_assign (s, NULL, GxB_PAIR_BOOL, tmasked, GrB_ALL, n,
            //    NULL) ;
            //LAGr_Vector_nvals (&ignore, s) ; // finish the work for assign
            // TODO why this is slow that the eWiseAdd above
            //LAGr_assign (s, tmasked, NULL, true, GrB_ALL, n, GrB_DESC_S) ;
            //GxB_print(s, 2);

            // t1 = LAGraph_toc(tic);
            // total_time3 += t1;

            // if nnz(tReq) == 0, no need to continue the rest of this loop
            GrB_Index tReq_nvals ;
            LAGr_Vector_nvals(&tReq_nvals, tReq);
            if (tReq_nvals == 0) { break; }

            // tless<tReq> = tReq .< t
            // TODO can tReq be a structural mask? see next line
            // If A(src,src) is not explicit 0 and A is all positive, then
            // yes
            // TODO currently assuming all edges weights > 0 (FIXME), so 
            // we're using a structural mask.  Explicit zeros in A would
            // require a valued mask.
            //printf("-------------------------------------------------\n");
            //GxB_print(tReq, 2);
            //GxB_print(t, 2);
            // LAGraph_tic (tic);
            // TODO: try clearing explicitly ...
            LAGr_Vector_clear (tless) ;
            LAGr_eWiseAdd (tless, tReq, NULL, GrB_LT_INT32, tReq, t,
                GrB_DESC_S  /* assumes all entries in A are > 0 */) ;
            // finish the work for assign
            //LAGr_Vector_nvals (&ignore, tless) ;

            // remove explicit zeros from tless so it can be used as a
            // structural mask
            GrB_Index tless_nvals ;
            LAGr_select (tless, NULL, NULL, GxB_NONZERO, tless, NULL, NULL);
            LAGr_Vector_nvals (&tless_nvals, tless) ;
            // t1 = LAGraph_toc(tic);
            // total_time4 += t1;
            if (tless_nvals == 0) { break ; }

            // update reachable node list/mask
            //LAGr_eWiseAdd (reach, NULL, NULL, GxB_PAIR_BOOL, reach, tless,
            //    NULL);
            LAGr_assign (reach, tless, NULL, true, GrB_ALL, n, GrB_DESC_S) ;

            // tmasked<tless> = select (i*delta <= tReq < (i+1)*delta)
            // since all entries of the 5 GAP graphs are known to be
            // positive, and the entries of tmasked are at least i*delta,
            // tReq = tmasked min.+ AL must be >= i*delta.
            // Therefore, there is no need to perform GxB_select with
            // GxB_GE_THUNK to find tmasked >= i*delta from tReq 
            // LAGraph_tic (tic);
            // better to clear vector before select (TODO: fix in GraphBLAS)
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
            // t1 = LAGraph_toc(tic);
            // total_time5 += t1;
            // GxB_print(tmasked, print_lvl);

            // t<tless> = tReq
            // GrB_apply is faster than GrB_eWiseAdd or GrB_assign here 
            // even when t is dense
            // LAGraph_tic (tic);
            // TODO use assign (but my apply should be fast too...)
            LAGr_apply (t, tless, NULL, GrB_IDENTITY_INT32, tReq,
                GrB_DESC_S) ;
            //LAGr_assign (t, tless, NULL, tReq, GrB_ALL, n, GrB_DESC_S) ;
            //LAGr_Vector_nvals (&ignore, t) ;
            // t1 = LAGraph_toc(tic);
            // total_time6 += t1;

            LAGr_Vector_nvals(&tmasked_nvals, tmasked);

        }

        // printf ("\n=============== next outer:\n") ;

        // tmasked<s> = t
        // LAGraph_tic (tic);
        LAGr_assign (tmasked, s, NULL, t, GrB_ALL, n, GrB_DESC_RS) ;
        // t1 = LAGraph_toc(tic);
        // total_time10 += t1;

        // tReq = AH'*tmasked
        // LAGraph_tic (tic);
        LAGr_vxm (tReq, NULL, NULL, GxB_MIN_PLUS_INT32, tmasked, AH, NULL) ;
        // t1 = LAGraph_toc(tic);
        // total_time2 += t1;

        // t = min(t, tReq)
        // When t is dense, it is best to get tless<tReq> = tReq .< t,
        // and use tless as mask to update t.
        //printf("----------------------------------------------------\n");
        //GxB_print(t, 2);
        // LAGraph_tic (tic);
        //------
        // best for sparse t:
        //LAGr_eWiseAdd(t, NULL, NULL, GrB_MIN_INT32, t, tReq, NULL);
        //------
        // best for dense t:  TODO use this instead when t is dense:
        LAGr_Vector_clear (tless) ;
        LAGr_eWiseAdd (tless, tReq, NULL, GrB_LT_INT32, tReq, t,
            GrB_DESC_S);
        LAGr_apply(t, tless, NULL, GrB_IDENTITY_INT32, tReq, NULL);
        //------
        // worse:
        //LAGr_eWiseAdd(t, tReq, NULL, GrB_MIN_INT32, t, tReq, NULL);
        // t1 = LAGraph_toc(tic);
        // total_time9 += t1;
        //GxB_print(t, 2);

        //------------------------------------------------------------------
        // find out how many left to be computed
        //------------------------------------------------------------------

        // LAGraph_tic (tic);
        // eWiseAdd is twice faster than assign TODO: why
        // update reachable node list/mask
        LAGr_assign (reach, tless, NULL, true, GrB_ALL, n, NULL) ;

        // remove previous buckets
        LAGr_assign (reach, s, NULL, false, GrB_ALL, n, GrB_DESC_S) ;
        LAGr_reduce (&remain, NULL, GxB_LOR_BOOL_MONOID, reach, NULL) ;
        // t1 = LAGraph_toc(tic);
        // total_time7 += t1;
        //GxB_print(reach, 2);

        LAGr_Vector_clear (s) ; // clear s for the next loop
        ++i;

    }

    /*
    printf("total time, select LT time, vxm time, update s, find tless,"
    "update tmasked, update t, select GE, update t time2, get t.*s time:\n"
    "%12.6g\n %12.6g\n %12.6g\n %12.6g\n %12.6g\n %12.6g\n"
    "%12.6g\n %12.6g\n %12.6g\n %12.6g \n",
    total_time, total_time1, total_time2, total_time3, total_time4,
    total_time5, total_time6, total_time7, total_time9, total_time10);
    */

    // result = t
    *path_length = t;
    t = NULL;

    // LAGraph_tic (tic) ;
    GrB_free (&lBound) ;
    GrB_free (&uBound) ;
    GrB_free (&Inf) ;
    GrB_free (&tmasked) ;
    GrB_free (&tReq) ;
    GrB_free (&tless) ;
    GrB_free (&s) ;
    GrB_free (&reach) ;
    GrB_free (&AL) ;
    GrB_free (&AH) ;
    // double t_free = LAGraph_toc (tic) ;

    // double total_time = LAGraph_toc(tic1);
    // printf("total time      %12.3f sec\n", total_time);
    // printf("init time       %12.3f sec, ratio %8.3f\n", t_pre, t_pre/total_time);
    // printf("select LT time  %12.3f sec, ratio %8.3f\n", total_time1, total_time1/total_time);
    // printf("vxm time        %12.3f sec, ratio %8.3f\n", total_time2, total_time2/total_time);
    // printf("update s time   %12.3f sec, ratio %8.3f\n", total_time3, total_time3/total_time);
    // printf("find tless time %12.3f sec, ratio %8.3f\n", total_time4, total_time4/total_time);
    // printf("update tmasked  %12.3f sec, ratio %8.3f\n", total_time5, total_time5/total_time);
    // printf("update t time   %12.3f sec, ratio %8.3f\n", total_time6, total_time6/total_time);
    // printf("select GE time  %12.3f sec, ratio %8.3f\n", total_time7, total_time7/total_time);
    // printf("update t time2  %12.3f sec, ratio %8.3f\n", total_time9, total_time9/total_time);
    // printf("tmasked<s>=t    %12.3f sec, ratio %8.3f\n", total_time10, total_time10/total_time);
// 
    // printf("free workspace  %12.3f sec, ratio %8.3f\n", t_free, t_free / total_time) ;
    return GrB_SUCCESS;
}
