
#define LG_FREE_ALL                             \
    printf ("done here: %d\n", __LINE__) ;      \
    printf ("msg: [%s]\n", msg) ;               \
    GrB_free (&centrality) ;                    \
    GrB_free (&A) ;                             \
    LAGraph_Delete (&G, msg) ;                  \

#include "LAGraphX.h"
#include "LG_internal.h"
#include <stdio.h>

double difference(GrB_Matrix bc, GrB_Matrix reference_bc)
{
    GrB_Matrix diff = NULL ;

    uint64_t n ;
    GrB_Matrix_nrows (&n, bc) ;

    // Compute diff = max(abs(reference_bc - bc))
    GrB_Matrix_new(&diff, GrB_FP64, n, n) ;
    GrB_eWiseAdd(diff, NULL, NULL, GrB_MINUS_FP64, reference_bc, bc, NULL) ;
    GrB_apply(diff, NULL, NULL, GrB_ABS_FP64, diff, NULL) ;

    double err = 1 ;
    GrB_reduce(&err, NULL, GrB_MAX_MONOID_FP64, diff, NULL) ;

    GrB_free(&diff) ;

    return err ;
} ;

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix centrality = NULL, A = NULL ;
    GrB_Info info ;

    // start GraphBLAS and LAGraph
    LAGRAPH_TRY (LAGraph_Init (msg)) ;

    //--------------------------------------------------------------------------
    // read in the graph via a Matrix Market file from stdin
    //--------------------------------------------------------------------------

    if (argc < 2)
    {
        printf("Usage: %s <matrix-market-file>\n", argv[0]);
        return (GrB_INVALID_VALUE) ;
    }

    FILE *f = fopen(argv[1], "r");
    if (f == NULL)
    {
        printf("Error: unable to open file %s\n", argv[1]);
        return (GrB_INVALID_VALUE) ;
    }

    double t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_MMRead (&A, f, msg)) ;
    fclose(f);
    uint64_t n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    LAGRAPH_TRY (LAGraph_DeleteSelfEdges (G, msg)) ;
    LAGRAPH_TRY (LAGraph_Cached_AT (G, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // compute edge betweenness centrality
    //--------------------------------------------------------------------------

    // LG_SET_BURBLE (true) ;

    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGr_EdgeBetweennessCentrality (&centrality, G, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGr_EdgeBetweennessCentrality: %g sec\n", t) ;

    // LG_SET_BURBLE (false) ;

    //--------------------------------------------------------------------------
    // check the results using LG_check_edgeBetweennessCentrality
    //--------------------------------------------------------------------------

    GrB_Matrix reference_centrality = NULL;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LG_check_edgeBetweennessCentrality(&reference_centrality, G, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LG_check_edgeBetweennessCentrality: %g sec\n", t) ;


    double err = difference(centrality, reference_centrality) ;
    printf ("Error between computed and reference centrality: %e\n", err) ;
    if (err < 1e-4)
    {
        printf ("Test passed.\n") ;
    }
    else
    {
        printf ("Test failure!\n") ;
    }

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    GrB_free (&centrality) ;
    GrB_free (&reference_centrality) ;
    LAGraph_Delete (&G, msg) ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
