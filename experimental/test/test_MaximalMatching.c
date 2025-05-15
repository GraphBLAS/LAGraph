//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_MaximalMatching
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Vidith Madhu, Texas A&M University

//------------------------------------------------------------------------------

/*
NOTE: Unlike the other tests, this does not use .mtx files, but rather generates the test
matrices using specified configurations and seeds with LAGraph_Random_Matrix

*/

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
    double matching_val ; // for unweighted matchings, the size of the matching. For weighted, sum of edge weights
    int matching_type ;   // 0: unweighted, 1: heavy, 2: light
    int is_exact ;        // whether or not matching_val is exactly computed for this test
    const char *name ;    // matrix file name
}
matrix_info ;

const matrix_info tests [ ] =
{
    // unweighted bipartite
    { 150, 0, 1, "random_unweighted_bipartite1.mtx" }, // 42
    { 150, 0, 1, "random_unweighted_bipartite2.mtx" }, // 69
    { 143, 0, 0, "random_unweighted_bipartite1.mtx" }, // repeat
    { 147, 0, 0, "random_unweighted_bipartite2.mtx" }, // repeat

    // unweighted general
    // { 25, 0, 1, 50, -1, 5.0 / 50, 31, "unweighted_general_1" } , // 31
    { 25, 0, 1, "random_unweighted_general1.mtx"},
    // { 100, 0, 1, 200, -1, 10.0 / 200, 101, "unweighted_general_2" }, // 101
    { 100, 0, 1, "random_unweighted_general2.mtx"},
    { 24, 0, 0, "random_unweighted_general1.mtx"},
    { 95, 0, 0, "random_unweighted_general2.mtx"},

    // weighted bipartite
    // answer, matching_type, is_exact, l_node, r_node, density, seed, name
    // { 3777422047635, 1, 0, 1000, 1000, 20.0 / 1000, 83, "weighted_bipartite_1" },
    { 775940425564, 1, 0, "random_weighted_bipartite1.mtx"}, // seed: 83, nodes: 500, spf: 8
    // { 9851292258178, 1, 0, 2500, 2500, 30.0 / 2500, 78, "weighted_bipartite_2" }
    { 417490248760, 1, 0, "random_weighted_bipartite2.mtx"}, // seed: 151, nodes: 300, spf: 5 
    { 181453589490, 2, 0, "random_weighted_bipartite1.mtx" }, // repeat
    { 133704435764, 2, 0, "random_weighted_bipartite2.mtx" }, // repeat
    // weighted general
    { 783685067769, 1, 0, "random_weighted_general1.mtx" }, // seed: 137, nodes: 500, spf: 8
    { 420609293186, 1, 0, "random_weighted_general2.mtx" }, // seed: 62, nodes: 300, spf: 5
    { 165090013148, 2, 0, "random_weighted_general1.mtx" }, // repeat
    { 128746478507, 2, 0, "random_weighted_general2.mtx" }, // repeat
    
    { 0, 0, 0, "" }
} ;

double thresholds [ ] = {
    0.85,   // unweighted matching, exact
    0.90,   // unweighted matching, naive
    0.80,   // weighted matching, naive, light
    0.90,   // weighted matching, naive, heavy
} ;

#define SEEDS_PER_TEST 10
#define LEN 512

char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

void test_MaximalMatching (void) 
{
#if LAGRAPH_SUITESPARSE
    OK (LAGraph_Init (msg)) ;
//  OK (LG_SET_BURBLE (true)) ;

    for (int k = 0 ; ; k++)
    {
        const char *aname = tests [k].name ;
        if (strlen (aname) == 0) break ;
        printf ("\n======================= %s:\n", aname) ;
        TEST_CASE (aname) ;

        // old code using files
        //--------------
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        TEST_MSG ("Filename %s is invalid", filename) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;
        //--------------

        TEST_CHECK (A != NULL) ;
        TEST_MSG ("Building of adjacency matrix failed") ;
        GxB_print (A, 1) ;

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

        OK (LAGraph_Cached_NSelfEdges (G, msg)) ;
        OK (LAGraph_Cached_AT (G, msg)) ;

//      if (G->nself_edges != 0)
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

        OK (LAGraph_Incidence_Matrix (&E, G, msg)) ;
        GrB_Index num_nodes ;
        GrB_Index num_edges ;
        OK (GrB_Matrix_nrows (&num_nodes, E)) ;
        OK (GrB_Matrix_ncols (&num_edges, E)) ;
        OK (GrB_Matrix_new (&E_t, GrB_FP64, num_edges, num_nodes)) ;
        OK (GrB_transpose (E_t, NULL, NULL, E, NULL)) ;
        // set to row major incase Incidence_Matrix gave col_major
        OK (GrB_set(E, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;

        // get weight vector
        OK (GrB_Vector_new (&weight, GrB_FP64, num_edges)) ;
        OK (GrB_reduce (weight, NULL, NULL, GrB_MAX_MONOID_FP64, E_t, NULL)) ;
        
        // used to check correctness of matching
        OK (GrB_Vector_new (&node_degree, GrB_UINT64, num_nodes)) ;

        OK (GrB_Vector_new (&hop_edges, GrB_BOOL, num_edges)) ;
        OK (GrB_Vector_new (&hop_nodes, GrB_BOOL, num_nodes)) ;

        double avg_acc = 0 ;
        size_t which_threshold ;
        uint64_t seed = 0 ;
        // run max matching
        for (int i = 0; i < SEEDS_PER_TEST; i++){
            // try random seeds
            OK (LAGraph_MaximalMatching (&matching, E, E_t, tests [k].matching_type, seed, msg)) ;
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
            double expected = tests [k].matching_val ;
            double matching_value = 0 ;

            if (tests [k].matching_type == 0) {
                // unweighted
                // we only care about the number of chosen edges
                uint64_t matching_val_int ;
                OK (GrB_Vector_nvals (&matching_val_int, matching)) ;
                matching_value = matching_val_int ;
                if (tests [k].is_exact) {
                    which_threshold = 0 ;
                } else {
                    which_threshold = 1 ;
                }
            } else {
                // weighted
                // sum the weights of the chosen edges.
                // eliminate entries in the weight that are not in the matching, then sum the remaining entries
                GrB_Vector use_weights = NULL ;

                OK (GrB_Vector_new (&use_weights, GrB_FP64, num_edges)) ;
                OK (GrB_eWiseMult (use_weights, NULL, NULL, GrB_TIMES_FP64, weight, matching, NULL)) ;
                OK (GrB_reduce (&matching_value, NULL, GrB_PLUS_MONOID_FP64, use_weights, NULL)) ;
                OK (GrB_free (&use_weights)) ;
                
                which_threshold = 3 ;
            }
            if (which_threshold == 0) {
                // exact matching must have strict upper bound
                printf ("matching_value %g expected %g\n", matching_value, expected) ;
                TEST_CHECK (matching_value <= expected) ;
            }
            double acc = matching_value / expected ;
            if (tests [k].matching_type == 2) {
                // flip it for light matchings
                acc = expected / matching_value ;
                which_threshold = 2 ;
            }
            avg_acc += acc ;
            OK (GrB_free (&matching)) ;
            seed += num_nodes ;
        }
        avg_acc /= SEEDS_PER_TEST ;
        ok = (avg_acc >= thresholds [which_threshold]) ;

        TEST_CHECK (ok) ;
        printf ("Value of produced matching has %.5f accuracy, passing threshold is %.5f\n for case (%d)\n", avg_acc, thresholds [which_threshold], k) ;

        OK (GrB_free (&A)) ;
        OK (GrB_free (&E)) ;
        OK (GrB_free (&E_t)) ;
        OK (GrB_free (&weight)) ;
        OK (GrB_free (&node_degree)) ;
        OK (GrB_free (&hop_edges)) ;
        OK (GrB_free (&hop_nodes)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }
    OK (LAGraph_Finalize (msg)) ;
#endif
}

void test_MaximalMatchingErrors (void)
{
#if LAGRAPH_SUITESPARSE
    OK (LAGraph_Init (msg)) ;
//  OK (LG_SET_BURBLE (true)) ;

    E = NULL ;
    matching = NULL ;

    OK (GrB_Matrix_new (&E, GrB_FP64, 1, 1)) ;

    // result pointer is null
    GrB_Info result = LAGraph_MaximalMatching (NULL, E, E, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // E matrix is null
    result = LAGraph_MaximalMatching (&matching, NULL, E, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // E_t matrix is null
    result = LAGraph_MaximalMatching (&matching, E, NULL, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    GrB_free (&E) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

TEST_LIST = {
    { "MaximalMatching", test_MaximalMatching },
    { "MaximalMatchingErrors", test_MaximalMatchingErrors },
    { NULL, NULL }
} ;
