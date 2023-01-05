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
    int matching_type ;
    bool is_exact ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    { 0, 0, "" }
} ;

double tolerances [ ] = {
    0.1,    // random matching, exact
    0,      // weighted matching, exact
    0.07,   // random matching, inexact
    0       // weighted matching, inexact
} ;

#define SEEDS_PER_TEST 10
#define LEN 512

char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

void test_MaximalMatching (void) 
{
    LAGraph_Init (msg) ;

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

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;

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
        // by definition, G->A must equal G->AT if G is undirected
        bool ok ;
        OK (LAGraph_Matrix_IsEqual (&ok, G->A, G->AT, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Input graph is not undirected") ;

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

        // run max matching
        for (int i = 0; i < SEEDS_PER_TEST; i++){
            // try random seeds
            uint64_t seed = rand() ;
            GrB_Vector matching = NULL ;
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
            ok = 0 ;
            size_t which_tolerance ;
            if (files [k].matching_type == 0) {
                // random
                // we only care about the number of chosen edges
                OK (GrB_Vector_nvals (&matching_value, matching)) ;
                if (files [k].is_exact) {
                    which_tolerance = 0 ;
                } else {
                    which_tolerance = 3 ;
                }
            } else {
                // weighted
                // sum the weights of the chosen edges.
                OK (GrB_reduce (&matching_value, matching, GrB_PLUS_MONOID_FP64, weight, GrB_DESC_S)) ;
                if (files [k].is_exact) {
                    which_tolerance = 1 ;
                } else {
                    which_tolerance = 4 ;
                }
            }
            double error = fabs(matching_value - expected) / expected ;
            ok = error <= tolerances [which_tolerance] ;
            TEST_CHECK (ok) ;
            TEST_MSG ("Value of produced matching has %.5f error for tolerance %.5f\n", error, tolerances [which_tolerance]) ;
        }    
    }
}

void test_MaximalMatchingErrors (void)
{

}

TEST_LIST = {
    {"MaximalMatching", test_MaximalMatching},
    {"MaximalMatchingErrors", test_MaximalMatchingErrors},
    {NULL, NULL}
} ;