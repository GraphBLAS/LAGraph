#define LG_FREE_WORK             \
    {                            \
        GrB_free(&T);            \
        GrB_free(&C);            \
        GrB_free(&W);            \
        GrB_free(&w_temp);       \
        \  
        GrB_free(&m);            \
        \   
        GrB_free(&m_index);      \
        GrB_free(&D);            \
        \                
        GrB_free(&E);            \
        GrB_free(&C_temp);       \
        \     
            GrB_free(&Identity); \
    \                                                                        
}

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
        GrB_free(C_f); \
    }

#include "LG_internal.h"

int LAGr_PeerPressureClustering(
    // output:
    GrB_Matrix *C_f, // Final clustering C_f[i][j] == 1 indicates vertex j is in cluster i
    // input:
    LAGraph_Graph G, // input graph
    char *msg)
{
    LG_CLEAR_MSG;

    GrB_Matrix T = NULL;

    // Cluster workspace matrix
    GrB_Matrix C = NULL;

    // Used for normalizing weights and assuring that vertices have equal votes
    GrB_Matrix W = NULL;

    GrB_Vector w_temp = NULL;

    GrB_Vector m = NULL;
    GrB_Vector m_index = NULL;
    GrB_Matrix D = NULL;
    GrB_Matrix E = NULL;

    GrB_Matrix C_temp = NULL;

    GrB_Matrix Identity = NULL;

    // Number of vertices in the graph
    GrB_Index n = 0;

    LG_ASSERT(C_f != NULL, GrB_NULL_POINTER);
    (*C_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix A = G->A;

    LG_ASSERT_MSG(G->out_degree != NULL, -106,
                  "G->out_degree must be defined");

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------
    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Matrix_new(&T, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&C_temp, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&D, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&E, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&Identity, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&C_f, GrB_BOOL, n, n));
    GRB_TRY(GrB_Vector_new(&w_temp, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m_index, GrB_INT64, n));

    // For now, assure that all vertices have equal weights
    LG_ASSERT_MSG(G->nself_edges == n, -106, "G->nself_edges must be equal to the number of nodes");

    // Normalize the weights by the out-degrees
    GRB_TRY(GrB_apply(w_temp, NULL, NULL, GrB_MINV_FP64, G->out_degree, NULL));
    GRB_TRY(GrB_Matrix_diag(&W, w_temp, 0));
    GRB_TRY(GrB_mxm(A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, W, A, GrB_DESC_R));

    // Initial cluster vector (each vertex is in its own cluster)
    GrB_Vector ones, trues;
    GrB_Index *ones_array, *ones_array_indices;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&ones_array, n, sizeof(GrB_UINT64), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&ones_array_indices, n, sizeof(GrB_Index), msg));

    GRB_TRY(GrB_Vector_new(&ones, GrB_UINT64, n));
    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Vector_new(&trues, GrB_BOOL, n));
    GRB_TRY(GrB_assign(trues, NULL, NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&Identity, trues, 0));
    GRB_TRY(GrB_Vector_free(&trues));

    GRB_TRY(GrB_Vector_extractTuples_UINT64(ones_array_indices, ones_array, &n, ones));

    GRB_TRY(GrB_Matrix_diag(&C, ones, 0)); // Initial cluster C_i
    GRB_TRY(GrB_Vector_free(&ones));

    GxB_print(W, GxB_COMPLETE);
    GxB_print(C, GxB_COMPLETE);

    GrB_Vector ones_fp;
    GRB_TRY(GrB_Vector_new(&ones_fp, GrB_FP64, n));
    GRB_TRY(GrB_assign(ones_fp, NULL, NULL, (double)1, GrB_ALL, n, NULL));

    //--------------------------------------------------------------------------
    // main algorithm logic
    //--------------------------------------------------------------------------
    while (true)
    {
        // T = C_i x A
        GRB_TRY(GrB_mxm(T, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C, A, GrB_DESC_R));

        // GRB_TRY (GrB_mxv(m, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP64, T, ones_fp, GrB_DESC_R)) ;

        // Maximum value vector
        GRB_TRY(GrB_mxv(m, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP64, T, ones_fp, GrB_DESC_RT0));

        // argmax code T. Davis SS User Guide p. 286
        // GRB_TRY(GrB_assign(m_index, NULL, NULL, 1, GrB_ALL, n, NULL));
        GRB_TRY(GrB_Matrix_diag(&D, m, 0));
        GRB_TRY(GrB_mxm(E, NULL, NULL, GxB_ANY_EQ_FP64, T, D, NULL));
        GRB_TRY(GrB_select(E, NULL, NULL, GrB_VALUENE_BOOL, E, 0, NULL));
        GRB_TRY(GrB_mxv(m_index, NULL, NULL, GxB_MIN_SECONDI_INT64, E, ones_fp, GrB_DESC_RT0));

        // m_index_values are ROW indices and m_index_indices are COLUMN indices
        GrB_Index *m_index_values, *m_index_indices;
        LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_values, n, sizeof(GrB_INT64), msg));
        LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_indices, n, sizeof(GrB_Index), msg));

        GRB_TRY(GrB_Vector_extractTuples_INT64(m_index_indices, m_index_values, &n, m_index));

        GRB_TRY(GrB_extract(C_temp, NULL, NULL, Identity, GrB_ALL, n, m_index_values, n, NULL));

        LAGraph_Free((void **)&m_index_values, NULL);
        LAGraph_Free((void **)&m_index_indices, NULL);

        bool res = NULL;
        LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&res, C, C_temp, msg));
        if (res)
        {
            GxB_print(T, GxB_COMPLETE);
            *C_f = C_temp; // Set output matrix
            break;
        }

        // Maybe unnecessary because of R descriptor?
        GRB_TRY(GrB_Matrix_dup(&C, C_temp));
        GRB_TRY(GrB_Matrix_clear(C_temp));
        GRB_TRY(GrB_Matrix_clear(T));
    }

    // GxB_print(T, GxB_COMPLETE);
    // GxB_print(m, GxB_COMPLETE);
    // GxB_print(D, GxB_COMPLETE);
    GxB_print(C_temp, GxB_COMPLETE);
    GxB_print(m_index, GxB_COMPLETE);

    GRB_TRY(GrB_Vector_free(&ones_fp));

    LG_FREE_ALL;

    return (GrB_SUCCESS);
}