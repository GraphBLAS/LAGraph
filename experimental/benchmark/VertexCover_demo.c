// Experimental demo for Vertex Cover

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"
#include "LG_Xtest.h"

// LG_FREE_ALL is required by LG_TRY
#undef LG_FREE_ALL
#define LG_FREE_ALL              \
    {                            \
        GrB_free(&VC);           \
        LAGraph_Delete(&G, msg); \
    }

int main(int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg[LAGRAPH_MSG_LEN]; // for error messages from LAGraph
    LAGraph_Graph G = NULL;
    GrB_Matrix A = NULL ;
    GrB_Vector VC = NULL; // Vertex Cover result vector
    GrB_Index n;

    // srand(time(NULL));
    // uint64_t seed = ((uint64_t)rand() << 32) | rand(); // seed for VC problem (Luby)

    // start GraphBLAS and LAGraph
    bool burble = false; // set true for diagnostic outputs
    demo_init(burble);
    LAGRAPH_TRY (LAGraph_Random_Init (msg)) ;

    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).

    double t = LAGraph_WallClockTime();
    char *matrix_name = (argc > 1) ? argv[1] : "stdin";
    LAGRAPH_TRY (readproblem (&G, NULL,
    true, true, true, NULL, false, argc, argv)) ;
    t = LAGraph_WallClockTime() - t;
    printf("Time to read the graph:      %g sec\n", t);

    printf("\n==========================The input graph matrix G:\n");
    // LG_TRY(LAGraph_Graph_Print(G, GxB_SHORT, stdout, msg));
    GxB_print(G->A, GxB_COMPLETE);

    // compute G->out_degree
    LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));

    // // compute G->in_degree, just to test it (not needed for any tests)
    // LAGRAPH_TRY(LAGraph_Cached_InDegree(G, msg));

    GRB_TRY(GrB_Matrix_nrows(&n, G->A));
    GRB_TRY(GrB_Vector_new(&VC, GrB_BOOL, n));

    //--------------------------------------------------------------------------
    // try the LAGraph_VertexCover algorithm
    //--------------------------------------------------------------------------

    // printf("self edges: %i", G->nself_edges) ;

    LAGRAPH_TRY( LAGraph_VertexCover(&VC, G, 1, NULL, msg) );
    // printf("Return code: %i", ret);

    LG_FREE_ALL;
    // LAGRAPH_TRY (LAGraph_Random_Finalize (msg)) ;
    LG_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS);
}