#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "LAGraph.h"

void test_PageRank(void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL, cmatlab = NULL, diff = NULL ;
    int niters = 0 ;

    // create the karate graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A
    OK (LAGraph_Property_RowDegree (G, msg)) ;

    // compute its pagerank
    OK (LAGr_PageRank (&centrality, &niters, G, 0.85, 1e-4, 100, msg)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    // compare with MATLAB: cmatlab = centrality (G, 'pagerank')
    float err = difference (centrality, karate_rank) ;
    printf ("\nkarate:   err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    LAGraph_Finalize (msg) ;
}
