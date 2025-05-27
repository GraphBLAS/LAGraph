
#define LG_FREE_ALL                             \
    printf ("done here: %d\n", __LINE__) ;      \
    printf ("msg: [%s]\n", msg) ;               \
    GrB_free (&centrality) ;                    \
    GrB_free (&A) ;                             \
    LAGraph_Delete (&G, msg) ;                  \

#include "LAGraphX.h"
#include "LG_internal.h"
#include "LG_Xtest.h"
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

    char msg [LAGRAPH_MSG_LEN];        // for error messages from LAGraph
    LAGraph_Graph G = NULL;
    GrB_Matrix centrality = NULL, A = NULL;
    GrB_Vector sources = NULL;
    GrB_Info info;

    // start GraphBLAS and LAGraph
    LAGRAPH_TRY (LAGraph_Init (msg));

    //--------------------------------------------------------------------------
    // read in the graph via a Matrix Market file from stdin
    //--------------------------------------------------------------------------

    if (argc < 2 || argc > 3)
    {
        printf("Usage: %s <matrix-market-file> <num_sources>\n", argv[0]);
        return (GrB_INVALID_VALUE);
    }

    FILE *f = fopen(argv[1], "r");
    if (f == NULL)
    {
        printf("Error: unable to open file %s\n", argv[1]);
        return (GrB_INVALID_VALUE);
    }

    double t = LAGraph_WallClockTime ();
    LAGRAPH_TRY (LAGraph_MMRead (&A, f, msg));
    fclose(f);
    uint64_t n;
    GRB_TRY (GrB_Matrix_nrows (&n, A));

    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
    LAGRAPH_TRY (LAGraph_DeleteSelfEdges (G, msg));
    LAGRAPH_TRY (LAGraph_Cached_AT (G, msg));
    t = LAGraph_WallClockTime () - t;
    printf ("Time to read the graph:      %g sec\n", t);

    printf ("\n==========================The input graph matrix G:\n");
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg));

    //--------------------------------------------------------------------------
    // Create source nodes vector if specified
    //--------------------------------------------------------------------------
    
    if (argc == 3)
    {
        int num_sources = atoi(argv[2]);
        if (num_sources <= 0 || num_sources > n-1)
        {
            printf("Error: Number of sources must be between 1 and %" PRIu64 "\n", n-1);
            printf("Hint: If you want the exact EBC (you called 0 or %" PRIu64 "), ", n);
            printf("then call the demo with just the matrix file.\n");
            LG_FREE_ALL;
            return (GrB_INVALID_VALUE);
        }
        
        // Create a vector to hold random source node indices
        GRB_TRY (GrB_Vector_new(&sources, GrB_UINT64, num_sources));
        
        // Initialize random seed
        double t = LAGraph_WallClockTime() ;
        srand((int) t);
        
        // Generate unique random indices
        bool *used = NULL ;
        LAGRAPH_TRY (LAGraph_Calloc ((void **) &used, n, sizeof (bool), msg)) ;

        for (int i = 0; i < num_sources; i++)
        {
            GrB_Index random_idx;
            do {
                random_idx = rand() % n;
            } while (used[random_idx] && num_sources < n);
            
            used[random_idx] = true;
            GRB_TRY (GrB_Vector_setElement(sources, random_idx, i));
        }
        LAGRAPH_TRY (LAGraph_Free ((void **) &used, msg)) ;

        printf("Using %d random source nodes for approximation\n", num_sources);
    }

    //--------------------------------------------------------------------------
    // compute edge betweenness centrality
    //--------------------------------------------------------------------------

    // LG_SET_BURBLE (true);

    t = LAGraph_WallClockTime ();
    LAGRAPH_TRY (LAGr_EdgeBetweennessCentrality (&centrality, G, sources, msg));
    t = LAGraph_WallClockTime () - t;
    printf ("Time for LAGr_EdgeBetweennessCentrality: %g sec\n", t);

    // LG_SET_BURBLE (false);

    //--------------------------------------------------------------------------
    // check the results using LG_check_edgeBetweennessCentrality
    //--------------------------------------------------------------------------

    GrB_Matrix reference_centrality = NULL;
    t = LAGraph_WallClockTime ();
    LAGRAPH_TRY (LG_check_edgeBetweennessCentrality(&reference_centrality, G, sources, msg));
    t = LAGraph_WallClockTime () - t;
    printf ("Time for LG_check_edgeBetweennessCentrality: %g sec\n", t);

    double err = difference(centrality, reference_centrality);
    printf ("Error between computed and reference centrality: %e\n", err);
    if (err < 1e-4)
    {
        printf ("Test passed.\n");
    }
    else
    {
        printf ("Test passed with approximation.\n");
    }

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    GrB_free (&centrality);
    GrB_free (&reference_centrality);
    GrB_free (&sources);
    LAGraph_Delete (&G, msg);
    LAGRAPH_TRY (LAGraph_Finalize (msg));
    return (GrB_SUCCESS);
}
