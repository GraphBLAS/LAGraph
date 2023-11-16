#include "LAGraph_demo.h"

#define LG_FREE_ALL                     \
{                                       \
    LAGraph_Delete (&G, msg) ;          \
    GrB_free (&A) ;                     \
}

int main (int argc, char**argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix C_f ;         // Clustering result vector

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    FILE *f ;

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        false, false, false, NULL, true, argc, argv)) ;

    // compute G->out_degree
    LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

    // compute G->in_degree, just to test it (not needed for any tests)
    LAGRAPH_TRY (LAGraph_Cached_InDegree (G, msg)) ;

    // compute G->nself_edges
    LAGRAPH_TRY (LAGraph_Cached_NSelfEdges (G, msg)) ;

    // GrB_Index n ;
    // GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    printf("Input Matrix:\n");
    GxB_print (G->A, GxB_COMPLETE);

    // Run Peer Pressure Clustering algorithm
    GRB_TRY(LAGr_PeerPressureClustering(&C_f, G, msg));

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
