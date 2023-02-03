#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

// #define dbg

const int DUMP = 1;

int main (int argc, char** argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix E = NULL ;
    GrB_Vector matching = NULL ;

    bool burble = false;
    demo_init (burble) ;

    LAGRAPH_TRY (LAGraph_Random_Init (msg)) ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, false, GrB_FP64, false, argc, argv)) ;
    
    GRB_TRY (LAGraph_A_to_E (&E, G, msg)) ;
    // LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;
    #ifdef dbg
        printf("printing E now: \n");
    #endif
    char name [LAGRAPH_MAX_NAME_LEN] ;
    LAGRAPH_TRY (LAGraph_Matrix_TypeName (name, E, msg)) ;
    #ifdef dbg
        printf("type of E is: %s\n", name) ;
        LAGRAPH_TRY (LAGraph_Matrix_Print (E, 4, stdout, msg)) ;
        printf("running max matching now...\n") ;
    #endif 
    LAGRAPH_TRY (LAGraph_MaximalMatching (&matching, E, 0, 5, msg)) ;
    LAGRAPH_TRY (LAGraph_Vector_Print (matching, DUMP ? LAGraph_COMPLETE : LAGraph_SHORT, stdout, msg)) ;
    if (DUMP) {
        // want to dump E to the file to verify correctness of matching
        // (we need E to know which nodes correspond to entries in the matching vector)
        LAGRAPH_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;
    }
    GRB_TRY (GrB_free (&matching)) ;
    return (GrB_SUCCESS) ;
}