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
// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University.
// Adapted from GraphBLAS Template Library (GBTL) by Scott McMillian.

// LAGraph_sssp computes the shortest path lengths from the specified source
// vertex to all other vertices in the graph.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_WORK                   \
{                                           \
    GrB_free(&AL);                          \
    GrB_free(&AH);                          \
    GrB_free(&t);                           \
    GrB_free(&tmasked);                     \
    GrB_free(&tReq);                        \
    GrB_free(&tBi);                         \
    GrB_free(&tnew);                        \
    GrB_free(&tcomp);                       \
    GrB_free(&tless);                       \
    GrB_free(&s);                           \
    GrB_free(&leq_delta);                   \
    GrB_free(&gt_delta);                    \
    GrB_free(&geq_idelta);                  \
    GrB_free(&select_in_range);             \
    GrB_free(&MIN_PLUS_FP64);               \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK;              \
    GrB_free(path_length);          \
}

double LAGRAPH_SSSP_THRESHOLD_VALUE = 0;
double LAGRAPH_SSSP_LB = 0;
double LAGRAPH_SSSP_UB = 0;

void geq_threshold(void *out, const void *in)
{
    const double in_dbl = *((const double*) in);
    bool out_bool = (in_dbl >= LAGRAPH_SSSP_THRESHOLD_VALUE);
    *((bool*) out) = out_bool;
}

void gt_threshold(void *out, const void *in)
{
    const double in_dbl = *((const double*) in);
    bool out_bool = (in_dbl > LAGRAPH_SSSP_THRESHOLD_VALUE);
    *((bool*) out) = out_bool;
}

void leq_threshold(void *out, const void *in)
{
    const double in_dbl = *((const double*) in);
    bool out_bool = (in_dbl <= LAGRAPH_SSSP_THRESHOLD_VALUE);
    *((bool*) out) = out_bool;
}

void in_range(void *out, const void *in)
{
    const double in_dbl = *((const double*) in);
    bool out_bool = (in_dbl <= LAGRAPH_SSSP_UB) && (in_dbl >= LAGRAPH_SSSP_LB);
    *((bool*) out) = out_bool;
}

GrB_Info LAGraph_sssp // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix graph,          // input graph, treated as if boolean in semiring
    GrB_Index source,          // source vertex from which to compute shortest paths
    double delta               // delta value for delta stepping
)
{
    (*path_length) = NULL;
    GrB_Index n = 0; // Number of nodes in the graph
    GrB_Index tcomp_nvals = 0;
    GrB_Index tmasked_nvals = 0;

    GrB_Matrix AL = NULL;
    GrB_Matrix AH = NULL;

    GrB_Vector t = NULL;
    GrB_Vector tmasked = NULL;
    GrB_Vector tReq = NULL;

    GrB_Vector tBi = NULL;
    GrB_Vector tnew = NULL;
    GrB_Vector tcomp = NULL;
    GrB_Vector tless = NULL;
    GrB_Vector s = NULL;

    GrB_UnaryOp leq_delta = NULL;
    GrB_UnaryOp gt_delta = NULL;
    GrB_UnaryOp geq_idelta = NULL;
    GrB_UnaryOp select_in_range = NULL;

    GrB_Semiring MIN_PLUS_FP64 = NULL;

    LAGr_UnaryOp_new(&leq_delta, &leq_threshold, GrB_BOOL, GrB_FP64);
    LAGr_UnaryOp_new(&gt_delta, &gt_threshold, GrB_BOOL, GrB_FP64);
    LAGr_UnaryOp_new(&geq_idelta, &geq_threshold, GrB_BOOL, GrB_FP64);
    LAGr_UnaryOp_new(&select_in_range, &in_range, GrB_BOOL, GrB_FP64);

    GrB_Semiring_new(&MIN_PLUS_FP64, GxB_MIN_FP64_MONOID, GrB_PLUS_FP64);
    
    LAGr_Matrix_nrows(&n, graph); // Get dimensions

    // Create the workspace vectors
    LAGr_Vector_new(&t, GrB_FP64, n);
    LAGr_Vector_new(&tmasked, GrB_FP64, n);
    LAGr_Vector_new(&tReq, GrB_FP64, n);

    LAGr_Vector_new(&tBi, GrB_BOOL, n);
    LAGr_Vector_new(&tnew, GrB_BOOL, n);
    LAGr_Vector_new(&tcomp, GrB_BOOL, n);
    LAGr_Vector_new(&tless, GrB_BOOL, n);
    LAGr_Vector_new(&s, GrB_BOOL, n);

    // t = infinity, t[src] = 0
    LAGr_Vector_setElement(t, 0, source);

    // AL = A .* (A <= delta)
    LAGr_Matrix_new(&AL, GrB_FP64, n, n);

    // The following requires a GraphBLAS 1.3-compliant implementation:
    //LAGr_apply(AL, GrB_NULL, GrB_NULL, GrB_LE_FP64, graph, delta, GrB_NULL);

    LAGRAPH_SSSP_THRESHOLD_VALUE = delta;
    LAGr_apply(AL, GrB_NULL, GrB_NULL, leq_delta, graph, GrB_NULL);

    LAGr_apply(AL, AL, GrB_NULL, GrB_IDENTITY_FP64, graph, LAGraph_desc_ooor);

    // AH = A .* (A > delta)
    LAGr_Matrix_new(&AH, GrB_FP64, n, n);

    // The following requires a GraphBLAS 1.3-compliant implementation:
    //LAGr_apply(AH, GrB_NULL, GrB_NULL, GrB_GT_FP64, graph, delta, GrB_NULL);

    LAGRAPH_SSSP_THRESHOLD_VALUE = delta;
    LAGr_apply(AH, GrB_NULL, GrB_NULL, gt_delta, graph, GrB_NULL);

    LAGr_apply(AH, AH, GrB_NULL, GrB_IDENTITY_FP64, graph, LAGraph_desc_ooor);

    GrB_Index i = 0;

    // t >= i*delta

    // The following requires a GraphBLAS 1.3-compliant implementation:
    //LAGr_apply(tcomp, GrB_NULL, GrB_NULL, GrB_GE_FP64, t, i*delta, GrB_NULL);

    LAGRAPH_SSSP_THRESHOLD_VALUE = i*delta;
    LAGr_apply(tcomp, GrB_NULL, GrB_NULL, geq_idelta, t, GrB_NULL);

    LAGr_apply(tcomp, tcomp, GrB_NULL, GrB_IDENTITY_BOOL, tcomp, GrB_NULL);

    LAGr_Vector_nvals(&tcomp_nvals, tcomp);

    // while (t >= i*delta) not empty
    while (tcomp_nvals > 0)
    {
        LAGr_Vector_clear(s);

        // tBi = t .* (i*delta <= t < (i+1)*delta)
        LAGRAPH_SSSP_LB = (i) * delta;
        LAGRAPH_SSSP_UB = (i+1.0) * delta;

        LAGr_apply(tBi, GrB_NULL, GrB_NULL, select_in_range, t, GrB_NULL);

        LAGr_apply(tBi, tBi, GrB_NULL, GrB_IDENTITY_FP64, tBi, LAGraph_desc_ooor);

        // tm<tBi> = t
        LAGr_apply(tmasked, tBi, GrB_NULL, GrB_IDENTITY_FP64, t, LAGraph_desc_ooor);

        LAGr_Vector_nvals(&tmasked_nvals, tmasked);

        while (tmasked_nvals > 0)
        {
            // tReq = AL' (min.+) (t .* tBi)
            LAGr_vxm(tReq, GrB_NULL, GrB_NULL, MIN_PLUS_FP64, tmasked, AL, GrB_NULL);

            // s = s + tBi
            LAGr_eWiseAdd(s, GrB_NULL, GrB_NULL, GrB_LOR, s, tBi, GrB_NULL);

            // tless<tReq> = tReq .< t
            LAGr_eWiseAdd(tless, tReq, GrB_NULL, GrB_LT_FP64, tReq, t, LAGraph_desc_ooor);

            // tBi<tless> = i*delta <= tReq < (i+1)*delta
            LAGr_apply(tBi, tless, GrB_NULL, select_in_range, tReq, LAGraph_desc_ooor);

            // t = min(t, tReq)
            LAGr_eWiseAdd(t, GrB_NULL, GrB_NULL, GrB_MIN_FP64, t, tReq, GrB_NULL);

            // tm<tBi> = t
            LAGr_apply(tmasked, tBi, GrB_NULL, GrB_IDENTITY_FP64, t, LAGraph_desc_ooor);

            LAGr_Vector_nvals(&tmasked_nvals, tmasked);
        }

        // (t .* s)
        LAGr_apply(tmasked, s, GrB_NULL, GrB_IDENTITY_FP64, t, LAGraph_desc_ooor);

        // tReq = AH'(t .* s)
        LAGr_vxm(tReq, GrB_NULL, GrB_NULL, MIN_PLUS_FP64, tmasked, AH, GrB_NULL);

        // t = min(t, tReq)
        LAGr_eWiseAdd(t, GrB_NULL, GrB_NULL, GrB_MIN_FP64, t, tReq, GrB_NULL);

        ++i;

        // t >= i*delta
        // The following line requires a GraphBLAS 1.3-compliant implementation
        //LAGr_apply(tcomp, GrB_NULL, GrB_NULL, GrB_GE_FP64, t, i*delta, GrB_NULL);
        LAGRAPH_SSSP_THRESHOLD_VALUE = i*delta;
        LAGr_apply(tcomp, GrB_NULL, GrB_NULL, geq_idelta, t, GrB_NULL);

        LAGr_apply(tcomp, tcomp, GrB_NULL, GrB_IDENTITY_BOOL, tcomp, LAGraph_desc_ooor);

        LAGr_Vector_nvals(&tcomp_nvals, tcomp);
    }

    // result = t
    LAGr_Vector_dup(path_length, t);

    LAGRAPH_FREE_WORK;
    return GrB_SUCCESS;
}
