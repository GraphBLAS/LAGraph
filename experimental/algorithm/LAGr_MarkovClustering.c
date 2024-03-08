#define LG_FREE_WORK           \
    {                          \
        GrB_free(&C);          \
        GrB_free(&C_temp);     \
        GrB_free(&vpc);        \
        GrB_free(&w);          \
        GrB_free(&D);          \
        GrB_free(&ones);       \
        GrB_free(&MSE);        \
        GrB_free(&argmax_v);   \
        GrB_free(&argmax_p);   \
        GrB_free(&zero_INT64); \
        GrB_free(&true_BOOL);  \
    }

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
        GrB_free(c_f); \
    }

#define DEBUG

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_MarkovClustering(
    // output:
    GrB_Vector *c_f, // output matrix C_f(i, j) == 1 means vertex j is in cluster i
    // input
    int e,                        // expansion coefficient
    int i,                        // inflation coefficient
    double pruning_threshold,     // threshold for pruning values
    double convergence_threshold, // MSE threshold for convergence
    int max_iter,                 // maximum iterations
    LAGraph_Graph G,              // input graph
    char *msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    GrB_Matrix C = NULL;      // Cluster workspace matrix
    GrB_Matrix C_temp = NULL; // The newly computed cluster matrix at the end of each loop
    GrB_Matrix CC = NULL;
    GrB_Vector vpc = NULL; // vertices per cluster

    GrB_Vector w = NULL; // weight vector to normalize C matrix

    GrB_Matrix D = NULL;    // Diagonal workspace matrix
    GrB_Vector ones = NULL; // Vector of all 1's, used primarily to create identity

    GrB_Matrix MSE = NULL; // Mean squared error between C and C_temp (between subsequent iterations)

    GrB_Vector argmax_v = NULL;
    GrB_Vector argmax_p = NULL;

    GrB_Scalar zero_INT64 = NULL;
    GrB_Scalar true_BOOL = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    GrB_Matrix A = G->A;       // Adjacency matrix of G
    GrB_Index n, nrows, ncols; // Dimension of A
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));
    GRB_TRY(GrB_Matrix_nrows(&ncols, A));

    LG_ASSERT(c_f != NULL, GrB_NULL_POINTER);
    (*c_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    LG_ASSERT_MSG(G->out_degree != NULL, -106,
                  "G->out_degree must be defined");

    LG_ASSERT_MSG(nrows == ncols, -1002,
                  "Input matrix must be square");

    n = nrows;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&C, GrB_FP32, n, n));
    GRB_TRY(GrB_Matrix_new(&CC, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&C_temp, GrB_FP32, n, n));
    GRB_TRY(GrB_Vector_new(&vpc, GrB_INT64, n));
    GRB_TRY(GrB_Matrix_new(&MSE, GrB_FP32, n, n));
    GRB_TRY(GrB_Matrix_new(&D, GrB_FP32, n, n));
    GRB_TRY(GrB_Vector_new(&w, GrB_FP32, n));
    GRB_TRY(GrB_Vector_new(&ones, GrB_FP32, n));
    GRB_TRY(GrB_Vector_new(&argmax_v, GrB_FP32, n));
    GRB_TRY(GrB_Vector_new(&argmax_p, GrB_INT64, n));
    GRB_TRY(GrB_Scalar_new(&zero_INT64, GrB_INT64));
    GRB_TRY(GrB_Scalar_new(&true_BOOL, GrB_BOOL));

    GRB_TRY(GrB_Scalar_setElement(zero_INT64, 0));
    GRB_TRY(GrB_Scalar_setElement(true_BOOL, 1));

    // Create identity
    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&D, ones, 0));

    // Add self-edge to each vertex
    if (G->nself_edges != n)
    {
        GRB_TRY(GrB_assign(A, A, NULL, D, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
        G->A = A;
        G->out_degree, G->in_degree = NULL;
        G->nself_edges = LAGRAPH_UNKNOWN;
        LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_InDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));
    }

    GRB_TRY(GrB_Matrix_dup(&C_temp, A));
    GRB_TRY(GrB_Matrix_dup(&C, C_temp));

    double mse = 0;
    GrB_Index iter = 0;
    GrB_Index nvals;

    while (true)
    {
        printf("Iteration %lu\n", iter);
        GxB_print(C_temp, GxB_SUMMARY);

        // Normalization step: Scale each column in C_temp to add up to 1
        // w = 1 ./ sum(A(:j))
        // D = diag(w)
        GRB_TRY(GrB_reduce(w, NULL, NULL, GrB_PLUS_MONOID_FP32, C_temp, GrB_DESC_RT0));
        GRB_TRY(GrB_apply(w, NULL, NULL, GrB_MINV_FP32, w, GrB_DESC_R));
        GRB_TRY(GrB_Matrix_diag(&D, w, 0));
        GRB_TRY(GrB_mxm(C_temp, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP32, C_temp, D, GrB_DESC_R));

        // Prune values less than some small threshold
        GRB_TRY(GrB_select(C_temp, NULL, NULL, GrB_VALUEGT_FP32, C_temp, pruning_threshold, NULL));

        // Compute mean squared error between subsequent iterations
        GRB_TRY(GxB_Matrix_eWiseUnion(MSE, NULL, NULL, GrB_MINUS_FP32, C_temp, zero_INT64, C, zero_INT64, NULL));
        GRB_TRY(GrB_eWiseMult(MSE, NULL, NULL, GrB_TIMES_FP32, MSE, MSE, NULL));
        GRB_TRY(GrB_reduce(&mse, NULL, GrB_PLUS_MONOID_FP32, MSE, NULL));
        GRB_TRY(GrB_Matrix_nvals(&nvals, C_temp));
        mse /= nvals;

#ifdef DEBUG
        printf("\tMSE at iteration %lu: %f\n", iter, mse);
        printf("\tCurrent size of cluster matrix (nvals): %lu\n", nvals);
#endif

        bool res = NULL;
        LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&res, C, C_temp, msg));
        if (res || iter > max_iter || mse < convergence_threshold)
        {
#ifdef DEBUG
            printf("\nTerminated after %lu iterations\n\n", iter);
#endif
            break;
        }

        // Set C to the previous iteration
        GRB_TRY(GrB_Matrix_dup(&C, C_temp));

        // Expansion step
        for (int i = 0; i < e - 1; i++)
        {
            GRB_TRY(GrB_mxm(C_temp, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP32, C_temp, C_temp, NULL));
        }

        // Inflation step
        GRB_TRY(GrB_Matrix_apply_BinaryOp2nd_FP32(C_temp, NULL, NULL, GxB_POW_FP32, C_temp, (double)i, NULL));

        iter++;
    }

    // argmax_v = max (C_temp) where argmax_v(j) = max (C_temp (:,j))
    GRB_TRY(GrB_mxv(argmax_v, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP32, C_temp, ones, GrB_DESC_T0));
    // D = daig (argmax_v)
    GRB_TRY(GrB_Matrix_diag(&D, argmax_v, 0));
    GRB_TRY(GrB_mxm(CC, NULL, NULL, GxB_ANY_EQ_FP32, C_temp, D, NULL));
    GRB_TRY(GrB_select(CC, NULL, NULL, GrB_VALUENE_BOOL, CC, 0, NULL));
    GRB_TRY(GrB_mxv(argmax_p, NULL, NULL, GxB_MIN_SECONDI_INT64, CC, ones, GrB_DESC_T0));

    // pi := array of argmax_p indices, px := array of argmax_p values
    GrB_Index *pi, *px = NULL;
    GrB_Index p_nvals;
    GRB_TRY(GrB_Vector_nvals(&p_nvals, argmax_p));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&pi, n, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&px, n, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(pi, px, &n, argmax_p));
    // Can reuse CC matrix, this will hold the final clustering result
    GRB_TRY(GrB_Matrix_clear(CC));
    GRB_TRY(GxB_Matrix_build_Scalar(CC, px, pi, true_BOOL, p_nvals));
    LAGraph_Free((void *)&pi, NULL);
    LAGraph_Free((void *)&px, NULL);

    GRB_TRY(GrB_reduce(vpc, NULL, NULL, GrB_PLUS_MONOID_INT64, CC, NULL));

#ifdef DEBUG
    printf("Vertices per cluster\n");
    GxB_print(vpc, GxB_SHORT);
#endif

    GrB_Index *CfI, *CfJ, ncc;

    GRB_TRY(GrB_Matrix_nvals(&ncc, CC));

    LAGRAPH_TRY(LAGraph_Malloc((void **)&CfI, ncc, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&CfJ, ncc, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Matrix_extractTuples_BOOL(CfI, CfJ, NULL, &ncc, CC));

    GrB_Vector c = NULL;
    GRB_TRY(GrB_Vector_new(&c, GrB_INT64, ncc));
    GRB_TRY(GrB_Vector_build(c, CfJ, CfI, ncc, NULL));

    GrB_Vector_wait(c, GrB_MATERIALIZE);

    LAGraph_Free((void **)&CfI, NULL);
    LAGraph_Free((void **)&CfJ, NULL);

    (*c_f) = c; // Set output matrix

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}
