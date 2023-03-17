#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>

GrB_Vector matching = NULL , weight = NULL, node_degree = NULL, hop_nodes = NULL, hop_edges = NULL ;
GrB_Matrix A = NULL, E = NULL, E_t = NULL ;
LAGraph_Graph G = NULL ;

typedef struct
{
    double matching_val ; // for unweighted matchings, the size of the set. For weighted, sum of edge weights
    int matching_type ;   // 0: unweighted, 1: heavy, 2: light
    bool is_exact ;       // whether or not matching_val is exactly computed for this test
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    // unweighted bipartite
    { 25, 0, 1, "random_bipartite_bool_1.mtx" },
    { 9, 0, 1, "random_bipartite_bool_2.mtx" },
    { 1550, 0, 1, "random_bipartite_bool_4.mtx" },
    { 2379, 0, 0, "random_bipartite_bool_5.mtx" },
    { 12290, 0, 0, "random_bipartite_bool_6.mtx" },
    // unweighted general
    { 4967, 0, 0, "random_general_bool_1.mtx" },
    { 8, 0, 1, "random_general_bool_2.mtx" },
    { 22, 0, 0, "random_general_bool_3.mtx" },
    { 2416, 0, 0, "random_general_bool_4.mtx" },
    { 499, 0, 1, "random_general_bool_5.mtx"},
    
    { 0, 0, 0, "" }
} ;

double thresholds [ ] = {
    0.85,   // random matching, exact
    0.93,   // random matching, naive
    0,      // weighted matching, naive, light
    0,      // weighted matching, naive, heavy
} ;

#define SEEDS_PER_TEST 10
#define LEN 512

char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

void test_MaximalMatching (void) 
{
    OK (LAGraph_Init (msg)) ;
    OK (LAGraph_Random_Init (msg)) ;

    for (int k = 0 ; ; k++)
    {
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        TEST_MSG ("Filename %s is invalid", filename) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;

        TEST_CHECK (A != NULL) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

        OK (LAGraph_Cached_NSelfEdges (G, msg)) ;
        OK (LAGraph_Cached_AT (G, msg)) ;

        if (G->nself_edges != 0)
        {
            // remove self-edges
            printf ("graph has %g self edges\n", (double) G->nself_edges) ;
            OK (LAGraph_DeleteSelfEdges (G, msg)) ;
            printf ("now has %g self edges\n", (double) G->nself_edges) ;
            TEST_CHECK (G->nself_edges == 0) ;
        }

        // check if G is undirected
        // by definition, G->A must equal G->AT iff G is undirected
        bool ok = 0;
        OK (LAGraph_Matrix_IsEqual (&ok, G->A, G->AT, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Input graph is not undirected") ;
        G->kind = LAGraph_ADJACENCY_UNDIRECTED ;

        OK (LAGraph_A_to_E (&E, G, msg)) ;
        GrB_Index num_nodes ;
        GrB_Index num_edges ;
        OK (GrB_Matrix_nrows (&num_nodes, E)) ;
        OK (GrB_Matrix_ncols (&num_edges, E)) ;
        OK (GrB_Matrix_new (&E_t, GrB_FP64, num_edges, num_nodes)) ;
        OK (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;

        // get weight vector
        OK (GrB_Vector_new (&weight, GrB_FP64, num_edges)) ;
        OK (GrB_reduce (weight, NULL, NULL, GrB_MAX_MONOID_FP64, E_t, NULL)) ;
        
        // used to check correctness of matching
        OK (GrB_Vector_new (&node_degree, GrB_UINT64, num_nodes)) ;

        OK (GrB_Vector_new (&hop_edges, GrB_BOOL, num_edges)) ;
        OK (GrB_Vector_new (&hop_nodes, GrB_BOOL, num_nodes)) ;

        printf("\n");
        double avg_slack = 0 ;
        size_t which_threshold ;
        // run max matching
        for (int i = 0; i < SEEDS_PER_TEST; i++){
            // try random seeds
            uint64_t seed = rand() ;
            OK (LAGraph_MaximalMatching (&matching, E, files [k].matching_type, seed, msg)) ;
            // check correctness
            OK (GrB_mxv (node_degree, NULL, NULL, LAGraph_plus_one_uint64, E, matching, NULL)) ;
            GrB_Index max_degree ;
            OK (GrB_reduce (&max_degree, NULL, GrB_MAX_MONOID_UINT64, node_degree, NULL)) ;
            TEST_CHECK (max_degree <= 1) ;
            TEST_MSG ("Matching is invalid") ;
            // check maximality
            // not maximal <-> there is an edge where neither node is matched
            // we can do a 1-hop from the matched edges and see if every edge in the graph is covered
            OK (GrB_mxv (hop_nodes, NULL, NULL, LAGraph_any_one_bool, E, matching, NULL)) ;
            OK (GrB_mxv (hop_edges, NULL, NULL, LAGraph_any_one_bool, E_t, hop_nodes, NULL)) ;
            GrB_Index hop_edges_nvals ;
            OK (GrB_Vector_nvals (&hop_edges_nvals, hop_edges)) ;
            TEST_CHECK (hop_edges_nvals == num_edges) ;
            TEST_MSG ("Matching is not maximal") ;
            // check that the value of the matching is close enough
            double expected = files [k].matching_val ;
            double matching_value = 0 ;

            if (files [k].matching_type == 0) {
                // random
                // we only care about the number of chosen edges
                uint64_t matching_val_int ;
                OK (GrB_Vector_nvals (&matching_val_int, matching)) ;
                matching_value = matching_val_int ;
                if (files [k].is_exact) {
                    which_threshold = 0 ;
                } else {
                    which_threshold = 1 ;
                }
            } else {
                // weighted
                // sum the weights of the chosen edges.
                // eliminate entries in the weight that are not in the matching, then sum the remaining entries
                OK (GrB_eWiseMult (weight, NULL, NULL, GrB_TIMES_FP64, weight, matching, NULL)) ;
                OK (GrB_reduce (&matching_value, NULL, GrB_PLUS_MONOID_FP64, weight, NULL)) ;
                
                which_threshold = 2 ;
            }
            if (which_threshold == 0) {
                // exact matching must have strict upper bound
                TEST_CHECK (matching_value <= expected) ;
            }
            double slack = matching_value / expected ;
            if (files [k].matching_type == 2) {
                // flip it for light matchings
                slack = expected / matching_value ;
            }
            avg_slack += slack ;
            OK (GrB_free (&matching)) ;
        }
        avg_slack /= SEEDS_PER_TEST ;
        ok = (avg_slack >= thresholds [which_threshold]) ;

        TEST_CHECK (ok) ;
        printf ("Value of produced matching has %.5f slack, tolerance is %.5f\n for case (%d)\n", avg_slack, thresholds [which_threshold], k) ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&E)) ;
        OK (GrB_free (&E_t)) ;
        OK (GrB_free (&weight)) ;
        OK (GrB_free (&node_degree)) ;
        OK (GrB_free (&hop_edges)) ;
        OK (GrB_free (&hop_nodes)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }
}

void test_MaximalMatchingErrors (void)
{
    OK (LAGraph_Init (msg)) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    TEST_MSG ("Loading of adjacency matrix failed") ;

    E = NULL ;
    matching = NULL ;

    OK (GrB_Matrix_new (&E, GrB_FP64, 1, 1)) ;

    // result pointer is null
    int result = LAGraph_MaximalMatching (NULL, E, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // E matrix is null
    result = LAGraph_MaximalMatching (&matching, NULL, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;
}

TEST_LIST = {
    { "MaximalMatching", test_MaximalMatching },
    { "MaximalMatchingErrors", test_MaximalMatchingErrors },
    { NULL, NULL }
} ;