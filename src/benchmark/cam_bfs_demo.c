#include "LAGraph_demo.h"

#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, msg) ;      \
    GrB_free (&A) ;                 \
    GrB_free (&Abool) ;             \
    GrB_free (&parent) ;            \
    GrB_free (&level) ;             \
    GrB_free (&SourceNodes) ;       \
}

int main (int argc, char**argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector level = NULL ;
    GrB_Vector parent = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    uint64_t seed = 1 ;
    FILE *f ;

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, &SourceNodes,
        false, false, true, NULL, false, argc, argv)) ;

    // compute G->out_degree
    LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

    // compute G->in_degree, just to test it (not needed for any tests)
    LAGRAPH_TRY (LAGraph_Cached_InDegree (G, msg)) ;

    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;


    GxB_print (G->A, GxB_COMPLETE);

    GrB_Index ntrials ;
    GRB_TRY (GrB_Matrix_nrows (&ntrials, SourceNodes)) ;

    int64_t src = 0;
    LG_BreadthFirstSearch_cam(&parent, G, src);
    GrB_free (&parent) ;

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;

}