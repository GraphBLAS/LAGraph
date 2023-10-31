

#define LG_FREE_WORK                         \
    {                                        \
        GrB_free(&S);                        \
        GrB_free(&frontier);                 \
        GrB_free(&paths);                    \
        GrB_free(&bc_update);                \
        GrB_free(&bc_score);                 \
        GrB_free(&weights);                  \      
        GrB_free(&workspace_1);              \   
        GrB_free(&workspace_2);              \                          
        GrB_free(&workspace_3);              \                          
        GrB_free(&workspace_4);              \                                               
    }

#define LG_FREE_ALL           \
    {                         \
        LG_FREE_WORK;         \
        GrB_free(centrality); \
    }

#include "LG_internal.h"

//--------------------------------------------------------------------------------------
// LAGr_Betweenness: vertex betweenness-centrality using Vanilla Brandes (not Batch)
//--------------------------------------------------------------------------------------

int LAGr_Betweenness_Brandes(
    // output:
    GrB_Vector *centrality, // centrality(i): betweeness centrality of i
    // input:
    LAGraph_Graph G,         // input graph
    // const GrB_Index *source, // Source vertex
    char *msg)
{

    LG_CLEAR_MSG;

    // The search, keeps track of the depth at which each vertex is seen
    GrB_Matrix S = NULL;

    // The number of shortest paths to each vertex
    // e.g., if paths[2] = 3, there are three shortest paths that go to vertex 3
    GrB_Vector paths = NULL;

    // The frontier, the number of shortest paths to vertices at the current depth
    GrB_Vector frontier = NULL;

    // The weights for the BC updates
    GrB_Vector weights = NULL;

    // The BC score for each vertex
    GrB_Vector bc_score = NULL;

    // The BC update for each vertex
    GrB_Vector bc_update = NULL;

    GrB_Vector workspace_1 = NULL;
    GrB_Vector workspace_2 = NULL;
    GrB_Vector workspace_3 = NULL;
    GrB_Vector workspace_4 = NULL;


    GrB_Index n = 0; // Number of nodes in the graph

    LG_ASSERT(centrality != NULL, GrB_NULL_POINTER);
    (*centrality) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix A = G->A;
    GrB_Matrix AT;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        AT = A;
    }
    else
    {
        // A and A' differ
        AT = G->AT;
        LG_ASSERT_MSG(AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required");
    }

    // =========================================================================
    // === initializations =====================================================
    // =========================================================================

    GRB_TRY(GrB_Matrix_nrows(&n, A)); // Get number of vertices
    GRB_TRY(GrB_Vector_new(&paths, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&frontier, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&bc_score, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&bc_update, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&weights, GrB_FP64, n));
    GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n, n));

    int64_t frontier_size = 0;

    GrB_Index depth;
    GrB_Index root ; // Current root value or starting vertex for which to compute BC updates 
    for (root = 0; root < n; root++)
    {
        printf("\n\n\n-------------------------------------------");
        depth = 0;
        GRB_TRY(GrB_Matrix_clear(S));
        GRB_TRY(GrB_Vector_clear(paths));
        int64_t one = 1;
        GRB_TRY(GrB_Vector_setElement(paths, one, root));

        // frontier = A(r, :)
        GRB_TRY(GrB_Col_extract(frontier, NULL, NULL, A, GrB_ALL, n, root, GrB_DESC_T0));
        GxB_print(A, GxB_COMPLETE);

        GRB_TRY(GrB_Vector_nvals(&frontier_size, frontier));
        while (frontier_size > 0)
        {
            printf("r = % i, Depth = %i\n", root, depth);
            GxB_print(frontier, GxB_COMPLETE);
            GxB_print(paths, GxB_COMPLETE);

            depth++;
            GRB_TRY(GrB_eWiseAdd(paths, NULL, NULL, GrB_PLUS_INT64, paths, frontier, NULL));
            GRB_TRY(GrB_assign(S, NULL, NULL, frontier, depth, GrB_ALL, n, NULL));
            GRB_TRY(GrB_vxm(frontier, paths, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, frontier, A, GrB_DESC_RC));

            GRB_TRY(GrB_Vector_nvals(&frontier_size, frontier));
        }
        while (depth >= 2)
        {
            GRB_TRY(GrB_Vector_new(&workspace_1, GrB_BOOL, n));
            GRB_TRY(GrB_Col_extract(workspace_1, NULL, NULL, A, GrB_ALL, n, depth, GrB_DESC_T0)); // workspace_1 = S(d,:)
            GRB_TRY(GrB_apply(workspace_1, NULL, NULL, GrB_PLUS_FP64, 1.0, workspace_1, NULL));     // workspace_2 =  
        }
        // printf("Printing S Matrix after algorithm\n");
        // GxB_print(S, GxB_COMPLETE);
    }



    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}