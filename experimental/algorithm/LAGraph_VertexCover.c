// Cam Quilici, Texas A&M
// Experimentation

#define debug

#define LG_FREE_WORK      \
    {                     \
        GrB_free(&mis);   \
        GrB_free(&cover); \     
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_VertexCover(
    // outputs:
    GrB_Vector *vc,
    // inputs:
    LAGraph_Graph G,        // input graph
    uint64_t seed,          // random number seed
    GrB_Vector ignore_node, // if NULL, no nodes are ignored.  Otherwise
                            // ignore_node(i) = true if node i is to be
                            // ignored, and not treated as a candidate
                            // added to maximal independent set.
    char *msg)
{
    LG_CLEAR_MSG;
    GrB_Vector mis = NULL;
    GrB_Vector cover = NULL;
    GrB_Matrix A; // G->A, the adjacency matrix
    GrB_Index n;  // # of nodes

    LG_TRY(LAGraph_CheckGraph(G, msg));
    LG_ASSERT(vc != NULL, GrB_NULL_POINTER);

    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        (G->kind == LAGraph_ADJACENCY_DIRECTED &&
         G->is_symmetric_structure == LAGraph_TRUE))
    {
        // the structure of A is known to be symmetric
        A = G->A;
    }
    else
    {
        // A is not known to be symmetric
        LG_ASSERT_MSG(false, -105, "G->A must be symmetric");
    }

    LG_ASSERT_MSG(G->out_degree != NULL, -106,
                  "G->out_degree must be defined");
    LG_ASSERT_MSG(G->nself_edges == 0, -107, "G->nself_edges must be zero");

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    // GRB_TRY(GrB_Vector_new(&mis, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&cover, GrB_BOOL, n));

    // srand(time(NULL));
    // uint64_t seed = ((uint64_t)rand() << 32) | rand();
    
    LAGRAPH_TRY( LAGraph_MaximalIndependentSet(&mis, G, (int64_t)1, NULL, msg)) ;
    // printf("Return code of MIS: %i\n", ret);

    GrB_assign(cover, mis, NULL, (bool)true, GrB_ALL, n, GrB_DESC_SC); // Structural complement

#ifdef debug
    printf("Maximal Independent Set:\n");
    LAGRAPH_TRY(LAGraph_Vector_Print(mis, LAGraph_SHORT, stdout, msg));
    printf("Minimal Vertex Cover:\n");
    LAGRAPH_TRY(LAGraph_Vector_Print(cover, LAGraph_SHORT, stdout, msg));
#endif

    GRB_TRY(GrB_wait(cover, GrB_MATERIALIZE));
    (*vc) = cover;
    LG_FREE_WORK;
    return (GrB_SUCCESS);
}