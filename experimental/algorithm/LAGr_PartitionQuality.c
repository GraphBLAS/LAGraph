#define LG_FREE_WORK                  \
    {                                 \
        GrB_free(&trace);             \
        GrB_free(&k);                 \
        GrB_free(&CA);                \
        GrB_free(&AT);                \
        GrB_free(&ONE_BOOL);          \
    }

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_PartitionQuality(
    // Outputs
    double *cov,  // Coverage
    double *perf, // Performance
    // Inputs
    GrB_Vector c, // Cluster vector where c[i] = j means vertex i is in cluster j
    GrB_Matrix A, // Adjacency matrix
    char *msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    GrB_Vector trace = NULL;
    GrB_Vector k = NULL;
    GrB_Matrix C = NULL;
    GrB_Matrix CA = NULL;

    GrB_Scalar ONE_BOOL = NULL;

    GrB_Matrix AT = NULL;

    GrB_Index n, nedges;
    uint64_t n_intraEdges, n_interEdges, n_interNonEdges, sum_k2;

    // Delete self-edges, not relevant to these clustering metrics
    GRB_TRY(GrB_select(A, NULL, NULL, GrB_OFFDIAG, A, 0, NULL));

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_nvals(&nedges, A));

    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&CA, GrB_INT64, n, n));
    GRB_TRY(GrB_Vector_new(&trace, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&k, GrB_INT64, n));
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

    // Check if A is symmetric
    bool is_symmetric;
    GRB_TRY(GrB_Matrix_new(&AT, GrB_BOOL, n, n));
    GRB_TRY(GrB_transpose(AT, NULL, NULL, A, NULL));
    LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&is_symmetric, A, AT, msg));

    // k = sum(C) .^ 2
    GRB_TRY(GrB_reduce(k, NULL, NULL, GrB_PLUS_MONOID_INT64, C, NULL));
    GRB_TRY(GrB_Vector_apply_BinaryOp2nd_INT64(k, NULL, NULL, GxB_POW_INT64, k, 2, NULL));
    // sum_k2 = total number of possible intra-cluster edges
    GRB_TRY(GrB_reduce(&sum_k2, NULL, GrB_PLUS_MONOID_INT64, k, NULL));

    // Calculate actual number of intra-cluster edges
    GRB_TRY(GrB_mxm(CA, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C, A, NULL));
    GxB_print(CA, GxB_COMPLETE);
    GRB_TRY(GrB_mxm(CA, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, CA, C, GrB_DESC_RT1));
    GxB_print(CA, GxB_COMPLETE);
    GRB_TRY(GxB_Vector_diag(trace, CA, 0, NULL));

    GxB_print(trace, GxB_SHORT);

    GRB_TRY(GrB_reduce(&n_intraEdges, NULL, GrB_PLUS_MONOID_INT64, trace, NULL));

    double coverage, performance;
    // Undirected
    if (is_symmetric)
    {
        nedges /= 2;
        n_intraEdges /= 2;
        n_interEdges = nedges - n_intraEdges;
        n_interNonEdges = (n * (n - 1) / 2) - ((sum_k2 - n) / 2) - n_interEdges;
        performance = (n_intraEdges + n_interNonEdges) * 1.0 / (n * (n - 1) / 2);
    }
    // Directed
    else
    {
        n_interEdges = nedges - n_intraEdges;
        // All possible edges - intra cluster non-edges gives space of possible
        // inter-cluster edges. Then subtract the actual inter-cluster edges to get
        // the number of inter-cluster non-edges.
        n_interNonEdges = n * (n - 1) - (sum_k2 - n) - n_interEdges;
        performance = (n_intraEdges + n_interNonEdges) * 1.0 / (n * (n - 1));
    }

    coverage = n_intraEdges * 1.0 / nedges;
    
    (*cov) = coverage;
    (*perf) = performance;

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}