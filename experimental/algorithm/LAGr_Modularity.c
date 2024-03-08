#define LG_FREE_WORK           \
    {                          \
        GrB_free(&l);          \
        GrB_free(&vmask);      \
        GrB_free(&k_in);       \
        GrB_free(&k_out);      \
        GrB_free(&in_degree);  \
        GrB_free(&out_degree); \
        GrB_free(&d);          \
        GrB_free(&CA);         \
        GrB_free(&ONE_BOOL);   \
    }

#define LG_FREE_ALL   \
    {                 \
        LG_FREE_WORK; \
    }

#include "LG_internal.h"
#include <LAGraphX.h>


int LAGr_Modularity(
    // Outputs
    double *mod_handle, // Modularity
    // Inputs
    double gamma, // Resolution parameter
    GrB_Vector c, // Cluster vector where c[i] = j means vertex i is in cluster j
    GrB_Matrix A, // Adjacency matrix
    char *msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    GrB_Vector l, vmask, k_in, k_out, d, out_degree, in_degree = NULL;
    GrB_Matrix C, CA = NULL;
    GrB_Scalar ONE_BOOL = NULL;

    GRB_TRY(GrB_select(A, NULL, NULL, GrB_OFFDIAG, A, 0, NULL));

    GrB_Index n, nedges;
    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_nvals(&nedges, A));

    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&CA, GrB_INT64, n, n));
    GRB_TRY(GrB_Vector_new(&l, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&vmask, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&k_in, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&k_out, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&d, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&out_degree, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&in_degree, GrB_INT64, n));
    GRB_TRY(GrB_Scalar_new(&ONE_BOOL, GrB_BOOL));

    GRB_TRY(GrB_Scalar_setElement_BOOL(ONE_BOOL, (bool)1));

    // Convert the cluster vector to a boolean matrix C where
    // C[i, j] = 1 if and only if vertex j is in cluster i
    GrB_Index *cI, *cX;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&cI, n, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&cX, n, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(cI, cX, &n, c));
    GRB_TRY(GxB_Matrix_build_Scalar(C, cX, cI, ONE_BOOL, n));
    GrB_Matrix_wait(C, GrB_MATERIALIZE);
    LAGraph_Free((void **)&cI, NULL);
    LAGraph_Free((void **)&cX, NULL);

    // Calculate actual number of intra-cluster edges
    GRB_TRY(GrB_mxm(CA, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C, A, NULL));
    GRB_TRY(GrB_mxm(CA, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, CA, C, GrB_DESC_RT1));
    GRB_TRY(GxB_Vector_diag(l, CA, 0, GrB_DESC_R));

    // Calculate the combined degree for each cluster
    GRB_TRY(GrB_reduce(out_degree, NULL, NULL, GrB_PLUS_MONOID_INT64, A, NULL));
    GRB_TRY(GrB_reduce(in_degree, NULL, NULL, GrB_PLUS_MONOID_INT64, A, GrB_DESC_T0));
    GRB_TRY(GrB_mxv(k_out, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C, out_degree, NULL));
    GRB_TRY(GrB_mxv(k_in, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C, in_degree, NULL));

    // vmask (i) = 0 if cluster i is non-empty (has aby vertices)
    GRB_TRY(GrB_reduce(vmask, NULL, NULL, GrB_LOR_MONOID_BOOL, C, NULL));
    GRB_TRY(GrB_apply(vmask, vmask, NULL, GxB_LNOT_BOOL, vmask, NULL));

    // If any of the above vectors have fewer entries than nclusters, this means that
    // there are singleton clusters with one vertex/no out-degree/no in-degree. So we 
    // need to add explicit zeros wherever values are missing for further calculations.
    GrB_Index nclusters, nl, nk_out, nk_in;
    GRB_TRY(GrB_Vector_nvals(&nclusters, vmask));
    GRB_TRY(GrB_Vector_nvals(&nl, l));
    GRB_TRY(GrB_Vector_nvals(&nk_out, l));
    GRB_TRY(GrB_Vector_nvals(&nk_in, l));

    if (nclusters != nl)
    {
        GRB_TRY(GrB_assign(l, l, NULL, vmask, GrB_ALL, nclusters, GrB_DESC_SC));
    }
    if (nclusters != nk_out)
    {
        GRB_TRY(GrB_assign(k_out, k_out, NULL, vmask, GrB_ALL, nclusters, GrB_DESC_SC));
    }
    if (nclusters != nk_in)
    {
        GRB_TRY(GrB_assign(k_in, k_in, NULL, vmask, GrB_ALL, nclusters, GrB_DESC_SC));
    }

    // Extract actual values of l, k_out, and k_in for modularity calculations
    GrB_Index *lX, *k_outX, *k_inX;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&lX, nclusters, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&k_outX, nclusters, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&k_inX, nclusters, sizeof(GrB_Index), msg));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, lX, &nclusters, l));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, k_outX, &nclusters, k_out));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, k_inX, &nclusters, k_in));

    GrB_Index m, out_degree_sum, in_degree_sum, L_c;
    GRB_TRY(GrB_reduce(&out_degree_sum, NULL, GrB_PLUS_MONOID_INT64, out_degree, NULL));

    m = out_degree_sum;
    double norm = 1.0 / (m * m);

    // Compute modularity
    double mod = 0.0;
    for (int c = 0; c < nclusters; c++)
    {
        mod += (1.0 * lX[c] / nedges) - (gamma * ((k_outX[c] * k_inX[c]) * norm));
    }

    (*mod_handle) = mod;

    LAGraph_Free((void **)&lX, NULL);
    LAGraph_Free((void **)&k_outX, NULL);
    LAGraph_Free((void **)&k_inX, NULL);

    (*mod_handle) = mod;

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}