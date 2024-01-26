//------------------------------------------------------------------------------
// LAGr_PeerPressureClustering: Graph clustering using the peer pressure method
//------------------------------------------------------------------------------

// Cameron Quilici, Texas A&M University

#define LG_FREE_WORK                     \
    {                                    \
        GrB_free(&T);                    \
        GrB_free(&C);                    \
        GrB_free(&W);                    \
        GrB_free(&w_temp);               \
        GrB_free(&m);                    \
        GrB_free(&m_index);              \
        GrB_free(&D);                    \
        GrB_free(&E);                    \
        GrB_free(&Identity_B);           \
        GrB_free(&Identity_F);           \
        GrB_free(&verts_per_cluster);    \
        GrB_free(&last_vpc);             \
        GrB_free(&diff_vpc);             \
        GrB_free(&zero_INT64);           \
        LAGraph_Free((void *)&AI, NULL); \
        LAGraph_Free((void *)&AJ, NULL); \
        LAGraph_Free((void *)&AX, NULL); \
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
    bool sanitize,
    LAGraph_Graph G, // input graph
    char *msg)
{
    // GxB_set(GxB_PRINT_1BASED, true);
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

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

    // Count number of vertices in each cluster
    GrB_Vector verts_per_cluster = NULL;
    GrB_Vector last_vpc = NULL;
    GrB_Vector diff_vpc = NULL;
    GrB_Scalar zero_INT64 = NULL;

    GrB_Matrix A_san = NULL;
    // Arrays holding extracted tuples if the matrix needs to be copied
    GrB_Index *AI = NULL;
    GrB_Index *AJ = NULL;
    double *AX = NULL;

    double p = .9;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    GrB_Matrix A;

    LG_ASSERT(C_f != NULL, GrB_NULL_POINTER);
    (*C_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Index n, nz;
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));
    GRB_TRY(GrB_Matrix_nvals(&nz, G->A));

    GrB_Vector ones_fp;
    GrB_Vector ones_fp_nz;
    GRB_TRY(GrB_Vector_new(&ones_fp, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&ones_fp_nz, GrB_FP64, nz));
    GRB_TRY(GrB_assign(ones_fp, NULL, NULL, (double)1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(ones_fp_nz, NULL, NULL, (double)1, GrB_ALL, nz, NULL));

    LAGraph_Matrix_TypeName(&MATRIX_TYPE, G->A, msg);
    printf("%s\n", MATRIX_TYPE);
    if (strcmp(MATRIX_TYPE, "float") != 0)
    {
        LAGRAPH_TRY(LAGraph_Malloc((void **)&AI, nz, sizeof(GrB_Index), msg));
        LAGRAPH_TRY(LAGraph_Malloc((void **)&AJ, nz, sizeof(GrB_Index), msg));
        GRB_TRY(GrB_Matrix_extractTuples_BOOL(AI, AJ, NULL, &nz, G->A));

        // Extract values from the ones_vector
        LAGRAPH_TRY(LAGraph_Malloc((void **)&AX, nz, sizeof(double), msg));
        GRB_TRY(GrB_Vector_extractTuples_FP64(NULL, AX, &nz, ones_fp_nz));

        // Build the sanitized matrix
        GRB_TRY(GrB_Matrix_new(&A_san, GrB_FP64, n, n));
        GRB_TRY(GrB_Matrix_build_FP64(A_san, AI, AJ, AX, nz, GrB_PLUS_FP64));

        A = A_san;
        GxB_print(A, GxB_SHORT);
    }
    else
    {
        A = G->A;
    }

    LG_ASSERT_MSG(G->out_degree != NULL, -106,
                  "G->out_degree must be defined");

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&T, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&C_temp, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&D, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&E, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&Identity_B, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&Identity_F, GrB_FP64, n, n)); // [?] Is this necessary?
    GRB_TRY(GrB_Vector_new(&w_temp, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m_index, GrB_INT64, n));

    GRB_TRY(GrB_Vector_new(&verts_per_cluster, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&last_vpc, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&diff_vpc, GrB_INT64, n));

    GRB_TRY(GrB_Scalar_new(&zero_INT64, GrB_INT64));
    GRB_TRY(GrB_Scalar_setElement(zero_INT64, 0));

    // Used below
    GrB_Vector ones, trues;
    GrB_Index *ones_array, *ones_array_indices;
    // [?] Is there a way to combine this call with the creation of the identity matrix creation
    // since in C boolean 1 is true? This feels unnecessary
    GRB_TRY(GrB_Vector_new(&ones, GrB_UINT64, n));
    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));

    GRB_TRY(GrB_Matrix_diag(&Identity_F, ones, 0)); // Identity FP Matrix

    // For now, assure that all vertices have equal weights
    // LG_ASSERT_MSG(G->nself_edges == n, -106, "G->nself_edges must be equal to the number of nodes");
    if (G->nself_edges != n)
    {
        GRB_TRY(GrB_assign(A, A, NULL, Identity_F, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
        G->A = A;

        G->out_degree = NULL;
        G->nself_edges = LAGRAPH_UNKNOWN;

        LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));
        // printf("nself edges %d\n", G->nself_edges);
#ifdef DEBUG
        GxB_print(A, GxB_SHORT);
#endif
    }

    //--------------------------------------------------------------------------
    // assuring vertices have equal votes by normalizing weights via out-degrees
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_apply(w_temp, NULL, NULL, GrB_MINV_FP64, G->out_degree, NULL));
    GRB_TRY(GrB_Matrix_diag(&W, w_temp, 0));
    GRB_TRY(GrB_mxm(A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, W, A, GrB_DESC_R));

    // Initial cluster vector (each vertex is in its own cluster)
    // Also initialize verts_per_cluster to be all 1 since each cluster starts with 1 vertex

    LAGRAPH_TRY(LAGraph_Malloc((void **)&ones_array, n, sizeof(GrB_UINT64), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void **)&ones_array_indices, n, sizeof(GrB_Index), msg));

    GRB_TRY(GrB_assign(verts_per_cluster, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(last_vpc, NULL, NULL, 1, GrB_ALL, n, NULL));

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

    GrB_Index last_num_changed = n;
    double tt, t0;

    //--------------------------------------------------------------------------
    // main algorithm logic
    //--------------------------------------------------------------------------
    GrB_Index count = 0;
    while (true)
    {

        // GxB_print(A, GxB_COMPLETE);
        tt = LAGraph_WallClockTime();
        // printf("Iteration %d\n", count);
        count++;

        // Tally (vote) matrix T where T[i][j] = k means there are k votes from cluster i for vertex j
        // to be in cluster i
        // T = C_i x A

        t0 = LAGraph_WallClockTime();
        GRB_TRY(GrB_mxm(T, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C, A, GrB_DESC_R));
        t0 = LAGraph_WallClockTime() - t0;
        printf("\tTime T = C * A (size = %i)\n\t%f\n", n, t0);

        // GRB_TRY(GrB_select(T, NULL, NULL, GrB_VALUENE_FP64, T, 0.0, NULL));

        // Maximum value vector where m[k] = l means l is the maximum fp value in column
        // k of the matrix T. In other words, the vector m holds the maximum number of votes that each
        // vertex got.

        t0 = LAGraph_WallClockTime();
        GRB_TRY(GrB_mxv(m, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP64, T, ones_fp, GrB_DESC_RT0));
        t0 = LAGraph_WallClockTime() - t0;
        printf("Time m = T * ones_fp (size = %i)\n\t%f\n", n, t0);

        // Now w\te need to find *which* cluster(s) cast the highest votes, for this we need argmax code
        // argmax ( adapted to argmin) code T. Davis SS User Guide p. 286
        t0 = LAGraph_WallClockTime();
        GRB_TRY(GrB_Matrix_diag(&D, m, 0));
        GRB_TRY(GrB_mxm(E, NULL, NULL, GxB_ANY_EQ_FP64, T, D, NULL));
        GRB_TRY(GrB_select(E, NULL, NULL, GrB_VALUENE_BOOL, E, 0, NULL)); // E == G in the pseudocode

        // m_index holds the first (i.e., if two clusters c_1 and c_2 cast the same vote L for vertex v
        // to gain membership into c_1/c_2, c_1 wins since it comes first/has the minimum index) row
        // index of T which is equal to the respective value of m, for each column
        GRB_TRY(GrB_mxv(m_index, NULL, NULL, GxB_MIN_SECONDI_INT64, E, ones_fp, GrB_DESC_RT0));

        // m_index_values are ROW indices and m_index_indices are COLUMN indices
        // [?] More of a C question but do these NEED to be allocated on heap? A more efficient way to do this?
        GrB_Index *m_index_values, *m_index_indices;
        LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_values, n, sizeof(GrB_INT64), msg));
        LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_indices, n, sizeof(GrB_Index), msg));

        GRB_TRY(GrB_Vector_extractTuples_INT64(m_index_indices, m_index_values, &n, m_index));

        GRB_TRY(GrB_extract(C_temp, NULL, NULL, Identity_B, GrB_ALL, n, m_index_values, n, NULL));
        t0 = LAGraph_WallClockTime() - t0;
        printf("\tArgmax time (size = %i)\n\t%f\n", n, t0);

        //--------------------------------------------------------------------------
        // begin debugging/printing/misc info section
        //--------------------------------------------------------------------------

        // Adds up the total number of vertices in each cluster, i.e., the number of nonzero entries in each
        // row of the current iteration of the cluster matrix
        GRB_TRY(GrB_reduce(verts_per_cluster, NULL, NULL, GrB_PLUS_MONOID_INT64, C_temp, GrB_DESC_R));

        // GxB_print(verts_per_cluster, GxB_SHORT);
        // GxB_print(last_vpc, GxB_SHORT);
        // GxB_print(diff_vpc, GxB_SHORT);

        LAGRAPH_TRY(GxB_eWiseUnion(diff_vpc, NULL, NULL, GrB_MINUS_INT64, verts_per_cluster, zero_INT64, last_vpc, zero_INT64, NULL));
        GRB_TRY(GrB_select(diff_vpc, NULL, NULL, GrB_VALUENE_INT64, diff_vpc, 0, GrB_DESC_R));

        GrB_Index num_changed = NULL;
        GRB_TRY(GrB_Vector_nvals(&num_changed, diff_vpc));

        // Re-use w_temp and W workspaces
        GRB_TRY(GrB_assign(w_temp, NULL, NULL, 1, GrB_ALL, n, GrB_DESC_R));
        // GxB_print(verts_per_cluster, GxB_COMPLETE);
//         GxB_print(w_temp, GxB_COMPLETE);
        GRB_TRY(GrB_apply(w_temp, verts_per_cluster, NULL, GxB_POW_FP64, verts_per_cluster, p, GrB_DESC_S));
        // GxB_print(w_temp, GxB_COMPLETE);

        GRB_TRY(GxB_Matrix_diag(W, w_temp, 0, GrB_DESC_R));
        // GxB_print(W, GxB_COMPLETE);
        // GxB_print(A, GxB_COMPLETE);
        GRB_TRY(GrB_mxm(A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, W, A, NULL));
        // GxB_print(A, GxB_COMPLETE);

        LAGraph_Free((void **)&m_index_values, NULL);
        LAGraph_Free((void **)&m_index_indices, NULL);

#ifdef DEBUG
        printf("\n--------------------------------------------------\n"
               "Current Values at iteration %i\n"
               "--------------------------------------------------\n",
               count);

        double percent_updated = num_changed * 1.0 / n;

        printf("Number of clusters updated since last iteration: %i\n", num_changed);
        printf("%2.3f %% of all cluster assignments have been updated since last iteration\n", percent_updated * 100);
        GxB_print(C_temp, GxB_SHORT);
        GxB_print(verts_per_cluster, GxB_SHORT);
        // GxB_print(last_vpc, GxB_SHORT);
        // GxB_print(diff_vpc, GxB_SHORT);
        GxB_print(m_index, GxB_SHORT);
        GxB_print(T, GxB_SHORT);
        printf("--------------------------------------------------\n\n\n");
#endif

        // Move back up eventually
        GRB_TRY(GrB_assign(last_vpc, NULL, NULL, verts_per_cluster, GrB_ALL, n, NULL));
        // last_num_changed = num_changed;

        //--------------------------------------------------------------------------
        // end debugging/printing/misc info section
        //--------------------------------------------------------------------------

        // When no changes to the cluster matrix have been made, terminate
        bool res = NULL;
        LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&res, C, C_temp, msg));
        if (res || count > 200)
        {
            (*C_f) = C_temp; // Set output matrix
            break;
        }

        // [?] Maybe unnecessary because of R descriptor? These DO need to be cleared
        // at the end of each iteration though...
        // Below code could potentially be very slow
        GRB_TRY(GrB_Matrix_dup(&C, C_temp));
        GRB_TRY(GrB_Matrix_clear(C_temp));
        GRB_TRY(GrB_Matrix_clear(T));

        tt = LAGraph_WallClockTime() - tt;
        printf("Total time of iteration %i (size = %i)\n\t%f\n\n\n", count, n, tt);
    }

    printf("--------------------------------------------------\n"
           "Final Information\n"
           "--------------------------------------------------\n"
           "Final tally matrix T where T[i][j] = k means there are "
           "k votes from cluster i for vertex j to be in cluster i:\n");
    GxB_print(T, GxB_SHORT);
    printf("Final cluster matrix C_temp where C_temp[i][j] == 1 means "
           "vertex j is in cluster i:\n");
    GxB_print(C_temp, GxB_SHORT);
    printf("Number of vertices per cluster:\n");
    GxB_print(verts_per_cluster, GxB_SHORT);

    GRB_TRY(GrB_Vector_free(&ones_fp));

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}