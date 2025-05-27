#include "../../src/benchmark/LAGraph_demo.h"   // for readproblem
#include "LAGraphX.h"                           // for LAGraph_coloring_independent_set
#include "LG_Xtest.h"                           // for LG_check_coloring
#include "LG_internal.h"                        // ?

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&C) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // setup variables, startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    // results
    GrB_Vector C = NULL ;    
    int num_colors = 0;
    double alg_time = 0;

    // other variables
    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;   

    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;
    LAGRAPH_TRY (LG_Random_Init (msg)) ;

    //--------------------------------------------------------------------------
    // read in the graph (defined in LAGraph_demo.h)
    //--------------------------------------------------------------------------

    double t = LAGraph_WallClockTime ( ) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        true,       // make the graph undirected, if true
        true,       // remove self-edges, if true
        true,       // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    // printf ("\n==========================The input graph matrix G:\n") ;
    // LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // execute independent set coloring algorithm
    //--------------------------------------------------------------------------
    
    t = LAGraph_WallClockTime ( ) ;
    int status = (LAGraph_coloring_independent_set (&C, &num_colors, G, msg)) ;
    alg_time = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for IS Coloring:        %g sec\n", alg_time) ;

    //--------------------------------------------------------------------------
    // check the results
    //--------------------------------------------------------------------------

    bool isequal ;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LG_check_coloring(G, C, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // print the results
    //--------------------------------------------------------------------------   
    
    printf ("\n===============================Time for IS       %g sec", alg_time) ;
    printf ("\n===============================Number of colors: %d\n\n", num_colors) ;
    
    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LG_Random_Finalize (msg)) ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;

}
