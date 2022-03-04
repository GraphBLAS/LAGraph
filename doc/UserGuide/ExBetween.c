#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "LAGraph.h"

void test_bc (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;
    int niters = 0 ;

    // create the karate graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A

    // compute its betweenness centrality
    OK (LAGr_Betweenness (&centrality, G, karate_sources, 4, msg)) ;
    printf ("\nkarate bc:\n") ;
    OK (LAGraph_Delete (&G, msg)) ;

    // compare with GAP:
    float err = difference (centrality, karate_bc) ;
    printf ("karate:   err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    LAGraph_Finalize (msg) ;
}

