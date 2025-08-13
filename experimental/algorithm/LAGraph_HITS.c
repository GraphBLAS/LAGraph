//------------------------------------------------------------------------------
// LAGraph_HITS: Hyperlink-Induced Topic Search algorithm using GraphBLAS API
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Aurko Routh, Texas A&M University

// TODO: almost ready for src.  this is advanced; add a basic method and add to src
// TODO: rename to LAGr_HubsAndAuthorities?

//------------------------------------------------------------------------------

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&h_old) ;             \
    GrB_free (&a_old) ;             \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free(&h);                   \
    GrB_free(&a);                   \
}

#include "LAGraphX.h"
#include "LG_internal.h"

// Check graph with error messages if it's empty

int LAGr_HITS
(
    GrB_Vector *hubs,
    GrB_Vector *authorities,
    int *iters,
    const LAGraph_Graph G,
    float tol,
    int itermax,
    char *msg
)
{

    LG_CLEAR_MSG ;
    GrB_Vector h = NULL, a = NULL, h_old = NULL, a_old=NULL;
    LG_ASSERT (hubs != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (authorities != NULL, GrB_NULL_POINTER) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GrB_Matrix AT ;
    GrB_Vector out_degree = G->out_degree, in_degree = NULL ;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        AT = G->A ;
        in_degree = out_degree ;
    }
    else
    {
        // A and A' differ
        AT = G->AT ;
        LG_ASSERT_MSG (AT != NULL,
            LAGRAPH_NOT_CACHED, "G->AT is required") ;
        in_degree = G->in_degree ;
    }
         // Initializations
    GrB_Index n;
    (*hubs) = NULL;
    (*authorities) = NULL;
    GRB_TRY(GrB_Matrix_nrows(&n, AT))
    float rdiff = 100;

    GRB_TRY (GrB_Vector_new (&h_old, GrB_FP32, n)) ;
    GRB_TRY (GrB_Vector_new (&a_old, GrB_FP32, n)) ;
    GRB_TRY (GrB_Vector_new (&h, GrB_FP32, n));
    GRB_TRY (GrB_Vector_new (&a, GrB_FP32, n)) ;

    float defaultValue = 1.0;
    GRB_TRY(GrB_assign(a, NULL, NULL, defaultValue, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(h, NULL, NULL, defaultValue, GrB_ALL, n, NULL));

    bool flag = true ;
    if (in_degree != NULL && out_degree != NULL)
    {
        GrB_Index in_nonempty = 0, out_nonempty = 0 ;
        // Count the number of non-zero entries in the in_degree vector
        GRB_TRY (GrB_Vector_nvals(&in_nonempty, in_degree)) ;
        // Count the number of non-zero entries in the out_degree vector
        GRB_TRY (GrB_Vector_nvals(&out_nonempty, out_degree)) ;
        flag = (in_nonempty + out_nonempty) > n/16.0;
    }

    for((*iters) = 0; (*iters) < itermax && rdiff > tol; (*iters)++) {
        // Save old values of h and a
        GrB_Vector temp = h_old ; h_old = h ; h = temp ;
        temp = a_old ; a_old = a ; a = temp ;

        if(flag) {
            //a = 0
            GRB_TRY(GrB_assign(a, NULL, GrB_PLUS_FP32, 0.0, GrB_ALL, n, NULL));
            //h = 0
            GRB_TRY(GrB_assign(h, NULL, GrB_PLUS_FP32, 0.0, GrB_ALL, n, NULL));
            // a += AT . h
            GRB_TRY(GrB_mxv(a, NULL,NULL, LAGraph_plus_second_fp32, AT, h_old, NULL));
            // h += A . a
            GRB_TRY(GrB_mxv(h, NULL,NULL, LAGraph_plus_second_fp32, G->A, a_old, NULL));
        } else {
            // a = AT . h
            GRB_TRY(GrB_mxv(a, NULL,NULL, LAGraph_plus_second_fp32, AT, h_old, NULL));
            // h = A . a
            GRB_TRY(GrB_mxv(h, NULL,NULL, LAGraph_plus_second_fp32, G->A, a_old, NULL));
        }

        // float maxA;
        float sumA;
        GRB_TRY(GrB_reduce(&sumA, NULL, GrB_PLUS_MONOID_FP32, a, NULL)); // Calculate the sum of all elements in the vector
        GRB_TRY(GrB_assign(a, NULL, GrB_DIV_FP32, sumA, GrB_ALL, n, NULL)); // Divide all elements by the sum

        float sumH;
        GRB_TRY(GrB_reduce(&sumH, NULL, GrB_PLUS_MONOID_FP32, h, NULL)); // Calculate the sum of all elements in the vector
        GRB_TRY(GrB_assign(h, NULL, GrB_DIV_FP32, sumH, GrB_ALL, n, NULL)); // Divide all elements by the sum

        // Deal with tolerance

        // a_old -= a
        GRB_TRY (GrB_assign(a_old, NULL, GrB_MINUS_FP32, a, GrB_ALL, n, NULL));

        // a_old = abs(a_old)
        GRB_TRY(GrB_apply (a_old, NULL, NULL, GrB_ABS_FP32, a_old, NULL));

        // rdiff = sum(a_old)
        GRB_TRY(GrB_reduce (&rdiff, NULL, GrB_PLUS_MONOID_FP32, a_old, NULL));

        // h_old -= h
        GRB_TRY (GrB_assign (h_old, NULL, GrB_MINUS_FP32, h, GrB_ALL, n, NULL));

        // h_old = abs(h_old)
        GRB_TRY (GrB_apply (h_old, NULL, NULL, GrB_ABS_FP32, h_old, NULL));

        // rdiff += sum(h_old)
        GRB_TRY (GrB_reduce (&rdiff, GrB_PLUS_FP32, GrB_PLUS_MONOID_FP32, h_old, NULL)) ;

        // rdiff = rdiff/2
        rdiff /= 2;
    }

    // Normalize
    float sumA;
    GRB_TRY(GrB_reduce(&sumA, NULL, GrB_PLUS_MONOID_FP32, a, NULL)); // Calculate the sum of all elements in the vector
    GRB_TRY(GrB_assign(a, NULL, GrB_DIV_FP32, sumA, GrB_ALL, n, NULL)); // Divide all elements by the sum

    float sumH;
    GRB_TRY(GrB_reduce(&sumH, NULL, GrB_PLUS_MONOID_FP32, h, NULL)); // Calculate the sum of all elements in the vector
    GRB_TRY(GrB_assign(h, NULL, GrB_DIV_FP32, sumH, GrB_ALL, n, NULL)); // Divide all elements by the sum

    (*hubs) = h;
    (*authorities) = a;
    LG_FREE_WORK;
    return (GrB_SUCCESS);
}
