//------------------------------------------------------------------------------
// LAGr_PeerPressureClustering: Graph clustering using the peer pressure method
//------------------------------------------------------------------------------

// Cameron Quilici, Texas A&M University

#define LG_FREE_WORK           \
    {                          \
        GrB_free(&T);          \
        GrB_free(&C);          \
        GrB_free(&W);          \
        GrB_free(&w_temp);     \
        GrB_free(&m);          \
        GrB_free(&m_index);    \
        GrB_free(&D);          \
        GrB_free(&E);          \
        GrB_free(&C_temp);     \
        GrB_free(&Identity_B); \
        GrB_free(&Identity_F); \
    }

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
        GrB_free(C_f); \
    }

// #define DEBUG

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_PeerPressureClustering(
    // output:
    GrB_Matrix *C_f, // Final clustering C_f[i][j] == 1 indicates vertex j is in cluster i
    // input:
    LAGraph_Graph G, // input graph
    char *msg)
{
    // GxB_set (GxB_PRINT_1BASED, true);

    GrB_Matrix T = NULL;

    // Cluster workspace matrix
    GrB_Matrix C = NULL;

    // The newly computed cluster matrix at the end of each loop
    GrB_Matrix C_temp = NULL;

    // Used for normalizing weights and assuring that vertices have equal votes
    GrB_Matrix W = NULL;
    GrB_Vector w_temp = NULL;

    // Objects used for the argmax functionality
    GrB_Vector m = NULL;
    GrB_Vector m_index = NULL;
    GrB_Matrix D = NULL;
    GrB_Matrix E = NULL;

    // The identity matrix of type GrB_BOOL
    GrB_Matrix Identity_B = NULL;
    GrB_Matrix Identity_F = NULL;

    // Number of vertices in the graph
    GrB_Index n = 0;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

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
    GRB_TRY(GrB_Matrix_new(&Identity_B, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&Identity_F, GrB_FP64, n, n)); // [?] Is this necessary?
    // GRB_TRY(GrB_Matrix_new(C_f, GrB_BOOL, n, n));
    GRB_TRY(GrB_Vector_new(&w_temp, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m_index, GrB_INT64, n));

    // Used below
    GrB_Vector ones, trues;
    GrB_Index *ones_array, *ones_array_indices;
    // [?] Is there a way to combine this call with the creation of the identity matrix creation
    // since in C boolean 1 is true? This feels unnecessary
    GRB_TRY(GrB_Vector_new(&ones, GrB_UINT64, n));
    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Matrix_diag(&Identity_F, ones, 0)); // Identity FP Matrix

    // For now, assure that all vertices have equal weights
    // printf("nselfedges %d", G->nself_edges);
    // LG_ASSERT_MSG(G->nself_edges == n, -106, "G->nself_edges must be equal to the number of nodes");
    if (G->nself_edges != n)
    {
        printf("Ensuring each vertex has a self edge\n");
        GRB_TRY(GrB_assign(A, A, NULL, Identity_F, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
        G->A = A;

        G->out_degree = NULL;
        G->nself_edges = LAGRAPH_UNKNOWN;

        LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));
        LAGRAPH_TRY (LAGraph_Cached_NSelfEdges (G, msg)) ;
        // printf("nself edges %d\n", G->nself_edges);
#ifdef DEBUG
        GxB_print(A, GxB_COMPLETE);
#endif
    }

    //--------------------------------------------------------------------------
    // assuring vertices have equal votes by normalizing weights via out-degrees
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_apply(w_temp, NULL, NULL, GrB_MINV_FP64, G->out_degree, NULL));
    GRB_TRY(GrB_Matrix_diag(&W, w_temp, 0));
    GRB_TRY(GrB_mxm(A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, W, A, GrB_DESC_R));

    // Initial cluster vector (each vertex is in its own cluster)

    LAGRAPH_TRY(LAGraph_Malloc((void **)&ones_array, n, sizeof(GrB_UINT64), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&ones_array_indices, n, sizeof(GrB_Index), msg));

    GRB_TRY(GrB_Vector_extractTuples_UINT64(ones_array_indices, ones_array, &n, ones));
    GRB_TRY(GrB_Matrix_diag(&C, ones, 0)); // Initial cluster C_i (this is the same as Identity with UINT_64 instead)
    GRB_TRY(GrB_Vector_free(&ones));

    // Create identity matrix of type boolean for later
    // [?] Maybe could avoid duplications and simply use the original C matrix above???
    GRB_TRY(GrB_Vector_new(&trues, GrB_BOOL, n));
    GRB_TRY(GrB_assign(trues, NULL, NULL, true, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&Identity_B, trues, 0));
    GRB_TRY(GrB_Vector_free(&trues));

    // GxB_print(W, GxB_COMPLETE);
    // GxB_print(C, GxB_COMPLETE);

    GrB_Vector ones_fp;
    GRB_TRY(GrB_Vector_new(&ones_fp, GrB_FP64, n));
    GRB_TRY(GrB_assign(ones_fp, NULL, NULL, (double)1, GrB_ALL, n, NULL));

    //--------------------------------------------------------------------------
    // main algorithm logic
    //--------------------------------------------------------------------------
    GrB_Index count = 0;
    while (true && count <= 20)
    {
        printf("Iteration %d\n", count);
        count++;
        // T = C_i x A
        GRB_TRY(GrB_mxm(T, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C, A, GrB_DESC_R));

        // Maximum value vector where m[k] = l means l is the maximum fp value in column
        // k of the matrix T
        GRB_TRY(GrB_mxv(m, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP64, T, ones_fp, GrB_DESC_RT0));

        // argmax code T. Davis SS User Guide p. 286
        GRB_TRY(GrB_Matrix_diag(&D, m, 0));
        GRB_TRY(GrB_mxm(E, NULL, NULL, GxB_ANY_EQ_FP64, T, D, NULL));
        GRB_TRY(GrB_select(E, NULL, NULL, GrB_VALUENE_BOOL, E, 0, NULL)); // E == G in the pseudocode
        // m_index holds the first row index of T which is equal to the respective value of m,
        // for each column
        GRB_TRY(GrB_mxv(m_index, NULL, NULL, GxB_MIN_SECONDI_INT64, E, ones_fp, GrB_DESC_RT0));

        // m_index_values are ROW indices and m_index_indices are COLUMN indices
        // [?] More of a C question but do these NEED to be allocated on heap?
        GrB_Index *m_index_values, *m_index_indices;
        LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_values, n, sizeof(GrB_INT64), msg));
        LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_indices, n, sizeof(GrB_Index), msg));

        GRB_TRY(GrB_Vector_extractTuples_INT64(m_index_indices, m_index_values, &n, m_index));

        GRB_TRY(GrB_extract(C_temp, NULL, NULL, Identity_B, GrB_ALL, n, m_index_values, n, NULL));

        LAGraph_Free((void **)&m_index_values, NULL);
        LAGraph_Free((void **)&m_index_indices, NULL);

        bool res = NULL;
        LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&res, C, C_temp, msg));
        if (res)
        {
            GxB_print(T, GxB_SHORT);
            *C_f = C_temp; // Set output matrix
            break;
        }

        // Unpack in order to get array indices
        // Apply keep as a dense vector

#ifdef DEBUG
        printf("--------------------------------------------------\n"
               "Current Values\n"
               "--------------------------------------------------\n");
        GxB_print(C_temp, GxB_COMPLETE);
        GxB_print(m_index, GxB_COMPLETE);
#endif

        // [?] Maybe unnecessary because of R descriptor? These do need to be cleared
        // at the end of each iteration though...
        GRB_TRY(GrB_Matrix_dup(&C, C_temp));
        GRB_TRY(GrB_Matrix_clear(C_temp));
        GRB_TRY(GrB_Matrix_clear(T));
    }

    printf("Final tally matrix T where T[i][j] = k means there are "
           "k votes from cluster i for vertex j to be in cluster i:\n");
    GxB_print(T, GxB_SHORT);
    printf("Final cluster matrix C_temp where C_temp[i][j] == 1 means "
           "vertex j is in cluster i:\n");
    GxB_print(C_temp, GxB_SHORT);
    GxB_print(m_index, GxB_SHORT);

    GRB_TRY(GrB_Vector_free(&ones_fp));

    LG_FREE_ALL;

    return (GrB_SUCCESS);
}