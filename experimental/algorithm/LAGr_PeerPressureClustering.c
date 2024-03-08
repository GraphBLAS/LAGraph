//------------------------------------------------------------------------------
// LAGr_PeerPressureClustering: Graph clustering using the peer pressure method
//------------------------------------------------------------------------------

// Cameron Quilici, Texas A&M University

#define LG_FREE_WORK                  \
    {                                 \
        GrB_free(&T);                 \
        GrB_free(&C);                 \
        GrB_free(&ones);              \
        GrB_free(&W);                 \
        GrB_free(&w_temp);            \
        GrB_free(&out_degree);        \
        GrB_free(&m);                 \
        GrB_free(&m_index);           \
        GrB_free(&D);                 \
        GrB_free(&E);                 \
        GrB_free(&I);                 \
    }

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
        GrB_free(c_f); \
    }

#define DEBUG

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_PeerPressureClustering(
    // output:
    GrB_Vector *c_f,      // Final clustering where c_f[i] = j means vertex i is in cluster j
    // input:
    bool normalize,       // if true, normalize the input graph via out-degree
    bool make_undirected, // if true, make G undirected which generally leads to a coarser partitioning
    double thresh,        // Threshold for convergence (percent of vertices that changed clusters)
    int max_iter,         // Maximum number of iterations
    LAGraph_Graph G,      // input graph
    char *msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    GrB_Matrix A = NULL;
    GrB_Matrix T = NULL;      // Tally matrix
    GrB_Matrix C = NULL;      // Cluster workspace matrix
    GrB_Matrix C_temp = NULL; // Newly computed cluster matrix at the end of each loop
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

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_ASSERT(c_f != NULL, GrB_NULL_POINTER);
    (*c_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    LG_ASSERT_MSG(G->out_degree != NULL, LAGRAPH_NOT_CACHED,
                  "G->out_degree is required");
    LG_ASSERT_MSG(G->AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required");

    GrB_Matrix A2;
    if (make_undirected && (G->kind == LAGraph_ADJACENCY_DIRECTED ||
                            G->is_symmetric_structure == LAGraph_FALSE))
    {
        // A and A' differ so set A = A + A'
        GRB_TRY(GrB_eWiseAdd(G->A, NULL, NULL, GrB_FIRST_FP64, G->A, G->AT, NULL));
    }
    A2 = G->A;

    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, A2));
    // All types of input matrices get cast to type FP64 for this algorithm
    GRB_TRY(GrB_Matrix_new(&A, GrB_FP64, n, n));
    GRB_TRY(GrB_apply(A, NULL, NULL, GrB_IDENTITY_FP64, A2, NULL));

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&T, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&C_temp, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&CD, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&D, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&E, GrB_BOOL, n, n));
    GRB_TRY(GrB_Matrix_new(&I, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&w_temp, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&m_index, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&out_degree, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&ones, GrB_FP64, n));


    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&I, ones, 0)); // Identity matrix of all 1 (float, bool, int)

    // Ensure that all vertices have self-edges
    GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_ONEB_FP64, A, I, NULL));

    //--------------------------------------------------------------------------
    // assuring vertices have equal votes by normalizing weights via out-degrees
    //--------------------------------------------------------------------------

    if (normalize)
    {
        GRB_TRY(GrB_reduce(out_degree, NULL, NULL, GrB_PLUS_MONOID_INT64, A, NULL));
        GRB_TRY(GrB_apply(w_temp, NULL, NULL, GrB_MINV_FP64, out_degree, NULL));
        GRB_TRY(GrB_Matrix_diag(&W, w_temp, 0));
        GRB_TRY(GrB_mxm(A, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, W, A, GrB_DESC_R));
    }

    // Initial cluster vector (each vertex is in its own cluster)
    GRB_TRY(GrB_Matrix_diag(&C, ones, 0));

    GrB_Index last_num_changed = n;
    GrB_Index num_changed;

    GrB_Index *m_index_values;
    LAGRAPH_TRY(LAGraph_Malloc((void **)&m_index_values, n, sizeof(GrB_INT64), msg));

    //--------------------------------------------------------------------------
    // main algorithm logic
    //--------------------------------------------------------------------------

    GrB_Index iter = 0;
    while (true)
    {
        // Tally (vote) matrix T where T[i][j] = k means there are k votes from cluster i for vertex j
        // to be in cluster i
        // T = C_i x A
        GRB_TRY(GrB_mxm(T, NULL, NULL, GxB_PLUS_SECOND_FP64, C, A, GrB_DESC_R));

        // Maximum value vector where m[k] = l means l is the maximum fp value in column
        // k of the matrix T. In other words, the vector m holds the maximum number of votes that each
        // vertex got.
        GRB_TRY(GrB_vxm(m, NULL, NULL, GrB_MAX_SECOND_SEMIRING_FP64, ones, T, GrB_DESC_R));

        // Now we need to find *which* cluster(s) cast the highest votes, for this we need argmax code
        // taken from T. Davis SS User Guide p. 286
        GRB_TRY(GrB_Matrix_diag(&D, m, 0));
        GRB_TRY(GrB_mxm(E, NULL, NULL, GxB_ANY_EQ_FP64, T, D, NULL));
        GRB_TRY(GrB_select(E, NULL, NULL, GrB_VALUENE_BOOL, E, 0, NULL)); // E == G in the pseudocode

        // m_index holds the first (i.e., if two clusters c_1 and c_2 cast the same vote L for vertex v
        // to gain membership into c_1/c_2, c_1 wins since it comes first/has the minimum index) row
        // index of T which is equal to the respective value of m, for each column
        GRB_TRY(GrB_vxm(m_index, NULL, NULL, GxB_MIN_SECONDI_INT64, ones, E, GrB_DESC_R));

        // m_index_values are row indices
        GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, m_index_values, &n, m_index));
        GRB_TRY(GrB_extract(C_temp, NULL, NULL, I, GrB_ALL, n, m_index_values, n, GrB_DESC_R));

        GRB_TRY(GrB_eWiseMult(CD, NULL, NULL, GrB_ONEB_BOOL, C, C_temp, GrB_DESC_R));
        GRB_TRY(GrB_reduce(&num_changed, NULL, GrB_PLUS_MONOID_INT64, CD, NULL));
        num_changed = n - num_changed;
        double percent_updated = num_changed * 1.0 / n;


        // When no changes to the cluster matrix have been made, terminate
        bool res = NULL;
        LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&res, C, C_temp, msg));
        if (res || percent_updated < thresh || iter > max_iter)
        {
            GrB_Index *CfI, *CfJ;
            LAGRAPH_TRY(LAGraph_Malloc((void **)&CfI, n, sizeof(GrB_Index), msg));
            LAGRAPH_TRY(LAGraph_Malloc((void **)&CfJ, n, sizeof(GrB_Index), msg));
            GRB_TRY(GrB_Matrix_extractTuples_BOOL(CfI, CfJ, NULL, &n, C_temp));

            GrB_Vector c = NULL;
            GRB_TRY(GrB_Vector_new(&c, GrB_INT64, n));
            GRB_TRY(GrB_Vector_build(c, CfJ, CfI, n, NULL));

            GrB_Vector_wait(c, GrB_MATERIALIZE);

            LAGraph_Free((void **)&CfI, NULL);
            LAGraph_Free((void **)&CfJ, NULL);

            (*c_f) = c; // Set output matrix
            break;
        }

        GRB_TRY(GrB_Matrix_dup(&C, C_temp));

        iter++;
    }

    LAGraph_Free((void **)&m_index_values, NULL);

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}
