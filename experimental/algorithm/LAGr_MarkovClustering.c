//------------------------------------------------------------------------------
// LAGr_MarkovClustering: Graph clustering using the Markov cluster algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Cameron Quilici, Texas A&M University

//------------------------------------------------------------------------------

// TODO: ready to consider for src

#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&T_prev);                                                     \
        GrB_free(&T);                                                          \
        GrB_free(&w);                                                          \
        GrB_free(&D);                                                          \
        GrB_free(&CC);                                                         \
        GrB_free(&ones);                                                       \
        GrB_free(&MSE);                                                        \
        GrB_free(&argmax_v);                                                   \
        GrB_free(&argmax_p);                                                   \
        GrB_free(&zero_FP32);                                                  \
    }

#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
        GrB_free(c_f);                                                         \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_MarkovClustering(
    // output:
    GrB_Vector *c_f, // output cluster vector
    // input
    int e,                        // expansion coefficient
    int i,                        // inflation coefficient
    double pruning_threshold,     // threshold for pruning values
    double convergence_threshold, // MSE threshold for convergence
    int max_iter,                 // maximum iterations
    LAGraph_Graph G,              // input graph
    char *msg)
{
#if LAGRAPH_SUITESPARSE
    GrB_Matrix T_prev = NULL;   // previous iteration transfer matrix
    GrB_Matrix T = NULL;        // current iteration transfer matrix
    GrB_Matrix CC = NULL;

    GrB_Vector w = NULL; // weight vector to normalize T_* matrix

    GrB_Matrix D = NULL;    // diagonal workspace matrix
    GrB_Vector ones = NULL; // vector of all 1

    GrB_Matrix MSE = NULL; // Mean squared error between T_prev and T (between
                           // subsequent iterations)

    GrB_Vector argmax_v = NULL;
    GrB_Vector argmax_p = NULL;

    GrB_Scalar zero_FP32 = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_ASSERT(c_f != NULL, GrB_NULL_POINTER);
    (*c_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Index n, nrows, ncols;
    GRB_TRY(GrB_Matrix_nrows(&nrows, G->A));
    GRB_TRY(GrB_Matrix_nrows(&ncols, G->A));

    LG_ASSERT_MSG(nrows == ncols, LAGRAPH_INVALID_GRAPH,
        "Input matrix must be square");
    n = nrows;
    LG_ASSERT_MSG(e >= 2, GrB_INVALID_VALUE, "e must be >= 2");

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&CC, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&MSE, GrB_FP32, n, n));
    GRB_TRY(GrB_Vector_new(&w, GrB_FP32, n));
    GRB_TRY(GrB_Vector_new(&ones, GrB_FP32, n));
    GRB_TRY(GrB_Vector_new(&argmax_v, GrB_FP32, n));
    GRB_TRY(GrB_Vector_new(&argmax_p, GrB_INT64, n));
    GRB_TRY(GrB_Scalar_new(&zero_FP32, GrB_FP32));
    GRB_TRY(GrB_Scalar_setElement(zero_FP32, 0));

    // Create identity
    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&D, ones, 0));

    // All types of input matrices get cast to type FP32 for this algorithm
    // and ensure all vertices have a self-edge
    GRB_TRY(GrB_Matrix_new(&T, GrB_FP32, n, n));
    GRB_TRY(GrB_eWiseAdd(T, NULL, NULL, GrB_FIRST_FP32, G->A, D, NULL));

    GrB_Index iter = 0;
    GrB_Index nvals;

    while (true)
    {
        // Normalization step: Scale each column in T to add up to 1
        // w = 1 ./ sum(T(:j))
        // D = diag(w)
        GRB_TRY(GrB_reduce(w, NULL, NULL, GrB_PLUS_MONOID_FP32, T,
                           GrB_DESC_T0));
        GRB_TRY(GrB_apply(w, NULL, NULL, GrB_MINV_FP32, w, NULL));
        GrB_free(&D);
        GRB_TRY(GrB_Matrix_diag(&D, w, 0));
        GRB_TRY(GrB_mxm(T, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP32,
                        T, D, NULL));

        if (iter > 0)
        {
            // Compute mean squared error between subsequent iterations
            double mse = 0;
            GRB_TRY(GxB_Matrix_eWiseUnion(MSE, NULL, NULL, GrB_MINUS_FP32, T,
                                          zero_FP32, T_prev, zero_FP32, NULL));
            GRB_TRY(GrB_apply(MSE, NULL, NULL, GxB_POW_FP32, MSE, (double) 2, NULL));

            GRB_TRY(GrB_reduce(&mse, NULL, GrB_PLUS_MONOID_FP32, MSE, NULL));
            GRB_TRY(GrB_Matrix_nvals(&nvals, MSE));
            mse /= nvals;
            if (iter > max_iter || mse < convergence_threshold)
                break;
        }

        // Set T_prev to the previous iteration
        GrB_free(&T_prev);
        T_prev = T ;
        T = NULL ;

        // compute the next T matrix
        // first expansion step for i = 0
        // T = T_prev^2
        GRB_TRY (GrB_Matrix_new (&T, GrB_FP32, n, n)) ;
        GRB_TRY(GrB_mxm(T, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP32,
                T_prev, T_prev, NULL));

        // subsequent expansion steps
        for (int i = 1; i < e - 1; i++)
        {
            // T = T^2
            GRB_TRY(GrB_mxm(T, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP32,
                            T, T, NULL));
        }

        // Inflation step
        GRB_TRY(GrB_Matrix_apply_BinaryOp2nd_FP32(
            T, NULL, NULL, GxB_POW_FP32, T, (double)i, NULL));

        // Prune values less than some small threshold
        GRB_TRY(GrB_select(T, NULL, NULL, GrB_VALUEGT_FP32, T,
                           pruning_threshold, NULL));

        iter++;
    }

    // In order to interpret the final iteration of the transfer matrix
    // (T), let an attractor vertex be a vertex with at least one positive
    // value within their corresponding row. Each attractor is attracting the
    // vertices (columns) which have positive values within its row. To compute
    // the output cluster vector, compute the argmax across columns of the
    // steady-state T matrix. Then argmax_p(i) = k means vertex i is in
    // cluster k.

    // argmax_v = max (T) where argmax_v(j) = max (T (:,j))
    GRB_TRY(GrB_mxv(argmax_v, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP32, T,
                    ones, GrB_DESC_T0));
    if (D != NULL) GrB_free(&D);
    GRB_TRY(GrB_Matrix_diag(&D, argmax_v, 0));
    GRB_TRY(GrB_mxm(CC, NULL, NULL, GxB_ANY_EQ_FP32, T, D, NULL));
    GRB_TRY(GrB_select(CC, NULL, NULL, GrB_VALUENE_BOOL, CC, 0, NULL));
    GRB_TRY(GrB_mxv(argmax_p, NULL, NULL, GxB_MIN_SECONDI_INT64, CC, ones,
                    GrB_DESC_T0));

    // Sometimes (particularly, when the pruning threshold is high), some
    // columns in the steady-state T have no values, i.e., they are not
    // attracted to any vertex. In this case, fill in the missing values with
    // the index of the vertex (these vertices will be arbitrarily put in the
    // cluster of their index).
    GrB_Index p_nvals;
    GRB_TRY(GrB_Vector_nvals(&p_nvals, argmax_p));
    if (p_nvals < n)
    {
        // argmax_p <!argmax_p> = 0:n-1
        GRB_TRY (GrB_apply (argmax_p, argmax_p, NULL, GrB_ROWINDEX_INT64,
            ones, 0, GrB_DESC_SC)) ;
    }

    (*c_f) = argmax_p ; // Set output vector
    argmax_p = NULL ;

    LG_FREE_WORK;

    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
