#include "LAGraphX.h"
#include "LG_internal.h"
#include <stdio.h>

#define LAGRAPH_CATCH(info)                     \
{                                               \
    GrB_free (&centrality) ;                    \
    GrB_free (&A) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
    return (info) ;                             \
}

int main (int argc, char **argv)
{
    double difference(GrB_Matrix bc, double* gap_result, GrB_Index rows, GrB_Index cols) ;

    double difference(GrB_Matrix bc, double* gap_result, GrB_Index rows, GrB_Index cols)
    {
        // GrB_Matrix diff = NULL;
        GrB_Matrix diff = NULL, gap_bc = NULL;
        OK(GrB_Matrix_new(&gap_bc, GrB_FP64, rows, cols));

        // Populate gap_bc with values from gap_result
        for (GrB_Index i = 0; i < rows; i++) {
            for (GrB_Index j = 0; j < cols; j++) {
                // if (*(gap_result + i * cols + j) != 0) printf("    (%ld, %ld)    %g\n", i, j, *(gap_result + i * cols + j));
                OK(GrB_Matrix_setElement_FP64(gap_bc, *(gap_result + i * cols + j), i, j));
            }
        }

        // Compute diff = max(abs(gap_bc - bc))
        OK(GrB_Matrix_new(&diff, GrB_FP64, rows, cols));
        OK(GrB_eWiseAdd(diff, NULL, NULL, GrB_MINUS_FP64, gap_bc, bc, NULL));
        OK(GrB_apply(diff, NULL, NULL, GrB_ABS_FP64, diff, NULL));

        double err = 0;
        OK(GrB_reduce(&err, NULL, GrB_MAX_MONOID_FP64, diff, NULL));

        OK(GrB_free(&diff));
        OK(GrB_free(&gap_bc));

        return err;
    }

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix centrality = NULL, A = NULL ;

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
    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // compute edge betweenness centrality
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGr_EdgeBetweennessCentrality (&centrality, G, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGr_EdgeBetweennessCentrality: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results using LG_check_edgeBetweennessCentrality
    //--------------------------------------------------------------------------

    GrB_Matrix reference_centrality = NULL;
    LAGRAPH_TRY (LG_check_edgeBetweennessCentrality(&reference_centrality, G, msg)) ;

    double err = difference(centrality, reference_centrality, G->n, G->n) ;
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
