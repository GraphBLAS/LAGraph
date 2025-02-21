#include "../../src/benchmark/LAGraph_demo.h"   // for readproblem
#include "LG_Xtest.h"                           // for LG_check_coloring
#include "LG_internal.h"                        // ?

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&C) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
    free(used_colors);                          \
    free(Ap);                                   \
    free(Ai);                                   \
    free(Ax);                                   \
}

void print_progress(double progress) {
    int bar_width = 50;
    printf("\r[");
    int pos = (int)(bar_width * progress);
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %.2f%%", progress * 100);
    fflush(stdout);
}

int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Vector C = NULL ;
    GrB_Matrix dupe = NULL;
    int *used_colors = NULL;
    GrB_Index *Ap = NULL;
    GrB_Index *Ai = NULL;
    void *Ax = NULL;
    int num_colors = 0;

    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read graph, unpack, and setup for algorithm
    //--------------------------------------------------------------------------

    // reading graph
    double t = LAGraph_WallClockTime ( ) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        true,       // make the graph undirected, if true
        false,       // remove self-edges, if true
        true,       // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;   

    // printf ("\n==========================The input graph matrix G:\n") ;
    // LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    // duplicate graph
    GrB_Matrix_dup (&dupe, G->A) ;

    // unpacking graph    
    GrB_Index Ap_size, Ai_size, Ax_size;
    GRB_TRY(GxB_Matrix_unpack_CSC(G->A, &Ap, &Ai, &Ax, &Ap_size, &Ai_size, &Ax_size, NULL, NULL, NULL));
    GrB_Index old_Ap_size = Ap_size;
    Ap_size = Ap_size / sizeof(GrB_Index);

    // setup
    GrB_Type Int = ((Ap_size - 1) < UINT32_MAX) ? GrB_UINT32 : GrB_UINT64 ;
    GrB_Vector_new (&C, Int, Ap_size - 1);
    used_colors = (int *) calloc (Ap_size - 1, sizeof(int));
    GrB_Index Ap_index;
    GrB_Index Ai_index;
    GrB_Index Ai_index_start, Ai_index_end;

    //--------------------------------------------------------------------------
    // execute greedy algorithm
    //--------------------------------------------------------------------------
    
    t = LAGraph_WallClockTime () ;   
    
    int current_color, neighbor_color;
    for (Ap_index = 0; Ap_index < Ap_size - 1; Ap_index++) {
        
        // setup bounds for Ai for this node
        // printf("setting bounds for node\n");
        Ai_index_start = Ap[Ap_index];
        Ai_index_end = Ap[Ap_index + 1];

        // for neighbors
        for (Ai_index = Ai_index_start; Ai_index < Ai_index_end; Ai_index++) {
            
            // skip self-edges
            if (Ai[Ai_index] == Ap_index) continue;
            
            // mark used colors
            if (GrB_Vector_extractElement(&neighbor_color, C, Ai[Ai_index]) == GrB_SUCCESS) {
                used_colors[neighbor_color] = 1;
            }
        }
                
        // find first unused color
        // printf("finding color for node \n");
        for (current_color = 1; current_color < Ap_size; current_color++) {
            if (used_colors[current_color] == 0) break;
        }
        if (current_color == Ap_size - 1) {
            printf("Error: all colors used\n");
            LG_FREE_ALL;
            return (1);
        }
        if (current_color > num_colors) num_colors = current_color;

        // assign color
        // printf("assigning color to node \n");
        GrB_Vector_setElement(C, current_color, Ap_index);
        
        // reset used colors
        // printf("resetting used colors\n");
        memset(used_colors, 0, sizeof(int) * (Ap_size - 1));

        // print progress
        print_progress(((double)Ap_index + 1) / (Ap_size - 1));
    }

    printf("\n");
    
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for Greedy Coloring:    %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results
    //--------------------------------------------------------------------------

    // replace graph
    G->A = dupe ;

    // // print G->A
    // printf("\n==========================The input graph matrix G:\n") ;
    // LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg) ;

    // // print C
    // printf("\n==========================The coloring vector C:\n") ;
    // LAGraph_Vector_Print (C, LAGraph_SHORT, stdout, msg) ;

    bool isequal ;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LG_check_coloring(G, C, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // print the results
    //--------------------------------------------------------------------------

    printf ("\n===============================The coloring vector C:\n") ;
    LAGraph_Vector_Print (C, LAGraph_SHORT, stdout, msg) ;

    printf ("\n===============================Time for Greedy:  %g sec", t) ;
    printf ("\n===============================Number of colors: %d\n\n", num_colors) ;
    
    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;

}
