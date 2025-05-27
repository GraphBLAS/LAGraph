//------------------------------------------------------------------------------
// LAGr_PeerPressureClustering: Graph clustering using the peer pressure method
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
// TODO: define the output vector c_f that defines the cluster assignment

#define LG_FREE_WORK                                                           \
    {                                                                          \
        GrB_free(&A);                                                          \
        GrB_free(&S);                                                          \
        GrB_free(&T);                                                          \
        GrB_free(&C);                                                          \
        GrB_free(&C_temp);                                                     \
        GrB_free(&CD);                                                         \
        GrB_free(&W);                                                          \
        GrB_free(&w_temp);                                                     \
        GrB_free(&out_degree);                                                 \
        GrB_free(&m);                                                          \
        GrB_free(&m_index);                                                    \
        GrB_free(&D);                                                          \
        GrB_free(&E);                                                          \
        GrB_free(&I);                                                          \
        GrB_free(&ones);                                                       \
        LAGraph_Free((void **)&m_index_values, NULL);                          \
        LAGraph_Free((void **)&CfI, NULL);                                     \
        LAGraph_Free((void **)&CfJ, NULL);                                     \
    }

#define LG_FREE_ALL                                                            \
    {                                                                          \
        LG_FREE_WORK;                                                          \
        GrB_free(c_f);                                                         \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_PeerPressureClustering(
    // output:
    GrB_Vector *c_f, // output cluster vector
    // input:
    bool normalize,       // if true, normalize the input graph via out-degree
    bool make_undirected, // if true, make G undirected which generally leads to
                          // a coarser partitioning
    double thresh,        // Threshold for convergence (percent of vertices that
                          // changed clusters)
    int max_iter,         // Maximum number of iterations
    LAGraph_Graph G,      // input graph
    char *msg)
{
#if LAGRAPH_SUITESPARSE

    GrB_Matrix A = NULL;
    GrB_Matrix S = NULL;      // symmetrized matrix, if needed
    GrB_Matrix T = NULL;      // Tally matrix
    GrB_Matrix C = NULL;      // Cluster workspace matrix
    GrB_Matrix C_temp = NULL; // Subsequent iteration cluster matrix
    GrB_Matrix CD = NULL;

    // Workspaces for weights (normalization and scaling)
    GrB_Matrix W = NULL;
    GrB_Vector w_temp = NULL;
    GrB_Vector out_degree = NULL;

    // Objects used for the argmax functionality
    GrB_Vector m = NULL;
    GrB_Vector m_index = NULL;
    GrB_Matrix D = NULL;
    GrB_Matrix E = NULL;

    // Identity matrix
    GrB_Matrix I = NULL;
    GrB_Vector ones = NULL;

    GrB_Index *m_index_values = NULL;
    GrB_Index *CfI = NULL, *CfJ = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_ASSERT(c_f != NULL, GrB_NULL_POINTER);
    (*c_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));

    GrB_Matrix A2;
    if (make_undirected && (G->kind == LAGraph_ADJACENCY_DIRECTED ||
                            G->is_symmetric_structure == LAGraph_FALSE))
    {
        // A and A' differ so set A = A + A'
        LG_ASSERT_MSG(G->AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required");
        GRB_TRY(GrB_Matrix_new(&S, GrB_FP64, n, n));
        GRB_TRY(GrB_eWiseAdd(S, NULL, NULL, GrB_ONEB_FP64, G->A, G->AT, NULL));
        A2 = S;
    }
    else
    {
        A2 = G->A;
    }

    // If the threshold is negative, set it to 0
    thresh = fmax(thresh, 0);

    // All types of input matrices get cast to type FP64 for this algorithm
    GRB_TRY(GrB_Matrix_new(&A, GrB_FP64, n, n));
    GRB_TRY(GrB_apply(A, NULL, NULL, GrB_IDENTITY_FP64, A2, NULL));

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&T, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&CD, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&E, GrB_BOOL, n, n));
    GRB_TRY(GrB_Vector_new(&m, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m_index, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&out_degree, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&ones, GrB_FP64, n));

    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));

    // Identity matrix of all 1 (cast throughout to float, bool, int)
    GRB_TRY(GrB_Matrix_diag(&I, ones, 0));

    // Ensure that all vertices have self-edges
    GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_ONEB_FP64, A, I, NULL));

    //--------------------------------------------------------------------------
    // assuring vertices have equal votes by normalizing weights via out-degrees
    //--------------------------------------------------------------------------

    if (normalize)
    {
        GRB_TRY(
            GrB_reduce(out_degree, NULL, NULL, GrB_PLUS_MONOID_INT64, A, NULL));
        GRB_TRY(GrB_Vector_new(&w_temp, GrB_FP64, n));
        GRB_TRY(GrB_apply(w_temp, NULL, NULL, GrB_MINV_FP64, out_degree, NULL));
        GRB_TRY(GrB_Matrix_diag(&W, w_temp, 0));
        GrB_free(&w_temp);
        GRB_TRY(
            GrB_mxm(A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, W, A, NULL));
    }

    // Initial cluster vector (each vertex is in its own cluster)
    GRB_TRY(GrB_Matrix_diag(&C, ones, 0));

    GrB_Index last_num_changed = n;
    GrB_Index num_changed;

    LG_TRY(LAGraph_Malloc((void **)&m_index_values, n, sizeof(GrB_Index), msg));

    //--------------------------------------------------------------------------
    // main algorithm logic
    //--------------------------------------------------------------------------

    GrB_Index iter = 0;
    while (true)
    {
        // Voting phase (T = A (plus,second) C)
        // T (i, j) = k <==> k votes for vertex j to be in cluster i
        GRB_TRY(GrB_mxm(T, NULL, NULL, GxB_PLUS_SECOND_FP64, C, A, NULL));

        // m = max (T (:, j))
        GRB_TRY(GrB_vxm(m, NULL, NULL, GrB_MAX_SECOND_SEMIRING_FP64, ones, T,
                        NULL));

        //------------------------------------------------------------------------
        // argmax across columns of T (T. Davis SS User Guide p. 286)
        //------------------------------------------------------------------------

        if (D != NULL)
            GrB_free(&D);
        GRB_TRY(GrB_Matrix_diag(&D, m, 0));
        GRB_TRY(GrB_mxm(E, NULL, NULL, GxB_ANY_EQ_FP64, T, D, NULL));
        // E = G in the pseudocode
        GRB_TRY(GrB_select(E, NULL, NULL, GrB_VALUENE_BOOL, E, 0, NULL));
        // Ties broken by minimum row index
        GRB_TRY(
            GrB_vxm(m_index, NULL, NULL, GxB_MIN_SECONDI_INT64, ones, E, NULL));
        // m_index_values(i) = argmax(T(:, i))
        GRB_TRY(
            GrB_Vector_extractTuples_INT64(NULL, (int64_t *) m_index_values, &n, m_index));
        GRB_TRY(GrB_Matrix_new(&C_temp, GrB_BOOL, n, n));
        GRB_TRY(GrB_extract(C_temp, NULL, NULL, I, GrB_ALL, n, m_index_values,
                            n, NULL));

        // Calculate change in cluster matrix between iterations
        GRB_TRY(GrB_eWiseMult(CD, NULL, NULL, GrB_ONEB_BOOL, C, C_temp, NULL));
        GRB_TRY(
            GrB_reduce(&num_changed, NULL, GrB_PLUS_MONOID_INT64, CD, NULL));
        num_changed = n - num_changed;
        double percent_updated = num_changed * 1.0 / n;

        // Terminate when cluster matrix reaches a steady-state
        if (percent_updated <= thresh || iter > max_iter)
        {
            // Convert cluster matrix to cluster vector
            LG_TRY(LAGraph_Malloc((void **)&CfI, n, sizeof(GrB_Index), msg));
            LG_TRY(LAGraph_Malloc((void **)&CfJ, n, sizeof(GrB_Index), msg));
            GRB_TRY(GrB_Matrix_extractTuples_BOOL(CfI, CfJ, NULL, &n, C_temp));

            GRB_TRY(GrB_Vector_new(c_f, GrB_INT64, n));
            GRB_TRY(GrB_Vector_build(*c_f, CfJ, CfI, n, NULL));

            GrB_Vector_wait(*c_f, GrB_MATERIALIZE);
            break;
        }

        GrB_free(&C);
        C = C_temp;
        C_temp = NULL;

        iter++;
    }

    LG_FREE_WORK;

    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
