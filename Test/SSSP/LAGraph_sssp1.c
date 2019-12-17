//------------------------------------------------------------------------------
// LAGraph_sssp: Single source shortest path with delta stepping
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

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
// LAGraph_sssp: Single source shortest path with delta stepping.
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

#include "sssp_test.h"

#define LAGRAPH_FREE_ALL                    \
{                                           \
    GrB_free(&AL);                          \
    GrB_free(&AH);                          \
    GrB_free(&t);                           \
    GrB_free(&tmasked);                     \
    GrB_free(&tReq);                        \
    GrB_free(&tBi);                         \
    GrB_free(&tnew);                        \
    GrB_free(&tless);                       \
    GrB_free(&s);                           \
}

// TODO assert the input matrix has type GrB_INT32
// TODO use INT32 instead of FP64
// TODO try different delta (for the GAP matrices)

GrB_Info LAGraph_sssp1         // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix graph,          // input graph, treated as if boolean in semiring (INT32)
    GrB_Index source,          // source vertex from which to compute shortest paths
    double delta               // delta value for delta stepping // TODO make int64 or int32
)
{
    GrB_Info info;

    (*path_length) = NULL;
    GrB_Index nrows, ncols, n = 0; // graph info
    GrB_Index tmasked_nvals = 0;

    GxB_Scalar lBound = NULL; // the threshold for GxB_select
    GxB_Scalar uBound = NULL; // the threshold for GxB_select

    GrB_Matrix AL = NULL;
    GrB_Matrix AH = NULL;

    GrB_Vector t = NULL;    // tentative shortest path length
    GrB_Vector tmasked = NULL;
    GrB_Vector tReq = NULL;

    GrB_Vector tBi = NULL;
    GrB_Vector tnew = NULL;
    GrB_Vector tless = NULL;
    GrB_Vector s = NULL;

    if (graph == NULL || path_length == NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    // Get dimensions
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, graph)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, graph)) ;

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

    LAGRAPH_OK (GxB_Scalar_new(&lBound, GrB_FP64));
    LAGRAPH_OK (GxB_Scalar_new(&uBound, GrB_FP64));
    LAGRAPH_OK (GxB_Scalar_setElement_FP64(lBound, delta));

    // Create the workspace vectors
    LAGr_Vector_new(&t, GrB_FP64, n);
    LAGr_Vector_new(&tmasked, GrB_FP64, n);
    LAGr_Vector_new(&tReq, GrB_FP64, n);

    LAGr_Vector_new(&tBi, GrB_BOOL, n);
    LAGr_Vector_new(&tnew, GrB_BOOL, n);
    LAGr_Vector_new(&tless, GrB_BOOL, n);
    LAGr_Vector_new(&s, GrB_BOOL, n);

    // t = infinity, t[src] = 0
    LAGr_Vector_setElement(t, 0, source);

    // AL = A .* (A <= delta) with lBound = delta
    LAGr_Matrix_new(&AL, GrB_FP64, n, n);
    LAGRAPH_OK (GxB_select(AL, GrB_NULL, GrB_NULL, GxB_LE_THUNK, graph, lBound,
        GrB_NULL));

    // AH = A .* (A > delta) with lBound = delta
    LAGr_Matrix_new(&AH, GrB_FP64, n, n);
    LAGRAPH_OK (GxB_select(AH, GrB_NULL, GrB_NULL, GxB_GT_THUNK, graph, lBound,
        GrB_NULL));

    GrB_Index i = 0;    // TODO use same INT32 or INT64

    // tmasked >= i*delta to find out how many left to be optimized
    LAGRAPH_OK (GxB_Scalar_setElement_FP64(lBound, i * delta));
    LAGRAPH_OK (GxB_select(tmasked, GrB_NULL, GrB_NULL, GxB_GE_THUNK, t, lBound,
        GrB_NULL));

    LAGr_Vector_nvals(&tmasked_nvals, tmasked);

    // while (t >= i*delta) not empty
    while (tmasked_nvals > 0)
    {
        LAGr_Vector_clear(s);

        // tmasked =  (tmasked < (i+1)*delta)
        LAGRAPH_OK (GxB_Scalar_setElement_FP64(uBound, (i+1.0) * delta));
        LAGRAPH_OK (GxB_select(tmasked, GrB_NULL, GrB_NULL, GxB_LT_THUNK,
            tmasked, uBound, GrB_NULL));

        // tBi = pattern of tmasked
        LAGr_apply(tBi, GrB_NULL, GrB_NULL, GxB_ONE_BOOL, tmasked,
            LAGraph_desc_ooor);

        LAGr_Vector_nvals(&tmasked_nvals, tmasked);

        while (tmasked_nvals > 0)
        {
            // tReq = AL' (min.+) (t .* tBi)
            LAGr_vxm(tReq, GrB_NULL, GrB_NULL, GxB_MIN_PLUS_FP64,
                tmasked, AL, GrB_NULL);

            // s = (s | tBi)
            LAGr_eWiseAdd(s, GrB_NULL, GrB_NULL, GrB_LOR, s, tBi,
                GrB_NULL);

            // tless<tReq> = tReq .< t
            LAGr_eWiseAdd(tless, tReq, GrB_NULL, GrB_LT_FP64, tReq,
                t, LAGraph_desc_ooor);

            // tmasked<tless> = i*delta <= tReq < (i+1)*delta
            // TODO try swapping the order; pick the fastest one
            LAGRAPH_OK (GxB_select(tmasked, tless, GrB_NULL, GxB_GE_THUNK,
                tReq,    lBound, GrB_NULL));
            LAGRAPH_OK (GxB_select(tmasked, tless, GrB_NULL, GxB_LT_THUNK,
                tmasked, uBound, GrB_NULL));

            // t<tless> = min(t, tReq)
            LAGr_eWiseAdd(t, tless, GrB_NULL, GrB_MIN_FP64, t,
                tReq, GrB_NULL);

            // tBi = pattern of tmasked
            LAGr_apply(tBi, GrB_NULL, GrB_NULL, GxB_ONE_BOOL,
                tmasked, LAGraph_desc_ooor);

            LAGr_Vector_nvals(&tmasked_nvals, tmasked);
        }

        // tmasked = (t .* s)
        LAGr_apply(tmasked, s, GrB_NULL, GrB_IDENTITY_FP64, t,
            LAGraph_desc_ooor);

        // tReq = AH'(t .* s)
        LAGr_vxm(tReq, GrB_NULL, GrB_NULL, GxB_MIN_PLUS_FP64,
            tmasked, AH, GrB_NULL);

        // t = min(t, tReq)
        LAGr_eWiseAdd(t, GrB_NULL, GrB_NULL, GrB_MIN_FP64, t, tReq,
            GrB_NULL);

        ++i;

        // t >= i*delta
        LAGRAPH_OK (GxB_Scalar_setElement_FP64(lBound, i * delta));
        LAGRAPH_OK (GxB_select(tmasked, GrB_NULL, GrB_NULL, GxB_GE_THUNK, t,
            lBound, GrB_NULL));

        LAGr_Vector_nvals(&tmasked_nvals, tmasked);
    }

    // result = t
    *path_length = t;
    t = NULL;

    LAGRAPH_FREE_ALL;
    return GrB_SUCCESS;
}
