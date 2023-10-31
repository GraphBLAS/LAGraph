#include "LAGraph_demo.h"

#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&centrality) ;        \
}

#define BRANDES 1

int main (int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Vector centrality = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        false, false, true, NULL, false, argc, argv)) ;
    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // set source nodes to be all vertices in the graph
    //--------------------------------------------------------------------------

    GrB_Index *vertex_list = malloc(n * sizeof(GrB_Index));
    if (!vertex_list) {
        printf("Failed to allocate memory for vertex_list\n");
        return (EXIT_FAILURE);
    }
    for (GrB_Index i = 0; i < n; i++) {
        vertex_list[i] = i;
    }

    //--------------------------------------------------------------------------
    // Compute betweenness centrality for all vertices
    //--------------------------------------------------------------------------

    double t_start = LAGraph_WallClockTime ( ) ;
    #if BRANDES == 1
    LAGRAPH_TRY (LAGr_Betweenness_Brandes (&centrality, G, msg)) ;
    #else
    LAGRAPH_TRY (LAGr_Betweenness (&centrality, G, vertex_list, n, msg)) ;
    #endif
    double t_end = LAGraph_WallClockTime ( ) ;

    printf ("BC time: %12.4f (sec)\n", t_end - t_start) ;
    fflush (stdout) ;

    GxB_print(centrality, GxB_COMPLETE);

    free(vertex_list);
    GrB_free (&centrality) ;

    LG_FREE_ALL;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
