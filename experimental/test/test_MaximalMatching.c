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

NOTE: Changes to LAGraph_Random may break these tests, since the LAGraph_Random implementation
used to build the test graphs may produce a different output from newer implementations
given the same seed.
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
    double matching_val ; // for unweighted matchings, the size of the set. For weighted, sum of edge weights
    int matching_type ;   // 0: unweighted, 1: heavy, 2: light
    int is_exact ;        // whether or not matching_val is exactly computed for this test
    GrB_Index n ;         // number of nodes in the graph. If bipartite, number of nodes in the left set
    GrB_Index m ;         // if not bipartite, should be -1. Otherwise, the number of nodes in the right set
    double density ;      // density of the matrix
    uint64_t seed ;       // seed used to generate the graph for this test
    const char *name ;
}
matrix_info ;

const matrix_info tests [ ] =
{
    // unweighted bipartite
    { 43, 0, 1, 50, 50, 5.0 / 50, 143, "random_bipartite_bool_1" },
    { 496, 0, 1, 500, 500, 3.0 / 500, 88, "random_bipartite_bool_2" },
    { 479, 0, 0, 500, 500, 10.0 / 500, 42, "random_bipartite_bool_3" },
    { 2483, 0, 0, 2500, 2500, 100.0 / 2500, 55, "random_bipartite_bool_4" },
    // unweighted general
    { 24, 0, 1, 50, -1, 5.0 / 50, 92, "random_general_bool_1" } ,
    { 100, 0, 1, 200, -1, 10.0 / 200, 112, "random_general_bool_2" },
    { 242, 0, 0, 500, -1, 10.0 / 500, 48, "random_general_bool_3" },
    { 1487, 0, 0, 3000, -1, 50.0 / 3000, 64, "random_general_bool_4" },
    // weighted bipartite
    { 3777422047635, 1, 0, 1000, 1000, 20.0 / 1000, 130, "random_bipartite_int_1" },
    { 9851292258178, 1, 0, 2500, 2500, 30.0 / 2500, 78, "random_bipartite_int_2" },
    { 372131180649, 2, 0, 1000, 1000, 20.0 / 1000, 24, "random_bipartite_int_3" },
    { 639851753175, 2, 0, 2500, 2500, 30.0 / 2500, 178, "random_bipartite_int_4" },
    // weighted general
    { 1847843295771, 1, 0, 1000, -1, 20.0 / 1000, 155, "random_general_int_1" },
    { 9991765577349, 1, 0, 5000, -1, 50.0 / 5000, 98, "random_general_int_2" },
    { 193597661237, 2, 0, 1000, -1, 20.0 / 1000, 44, "random_general_int_3" },
    { 520480326025, 2, 0, 5000, -1, 50.0 / 5000, 101, "random_general_int_4" },
    
    { 0, 0, 0, 0, 0, 0.0, 0, "" }
} ;

double thresholds [ ] = {
    0.85,   // random matching, exact
    0.90,   // random matching, naive
    0.80,   // weighted matching, naive, light
    0.90,   // weighted matching, naive, heavy
} ;

#define SEEDS_PER_TEST 10
#define LEN 512

char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

void test_MaximalMatching (void) 
{
    OK (LAGraph_Init (msg)) ;
//  GrB_set (GrB_GLOBAL, (int32_t) (true), GxB_BURBLE) ;
    OK (LAGraph_Random_Init (msg)) ;

    for (int k = 0 ; ; k++)
    {
        const char *aname = tests [k].name ;
        if (strlen (aname) == 0) break ;
        TEST_CASE (aname) ;

        // graph generation below
        //--------------
        if (tests [k].m != -1) {
            // we want a bipartite graph
            GrB_Index n = tests [k].n ;
            GrB_Index m = tests [k].m ;
            /*
            Bipartite graph generation works as follows: Create a random n x m matrix for the top right
            quadrant. The bottom left quadrant will be the transpose of this matrix. The other 2 quadrants
            will be empty.
            This exactly matches the process that the custom tests use
            */
            GrB_Matrix A_tr = NULL ; // top-right quadrant
            GrB_Index A_tr_nvals ;   // number of entries in top-right quadrant
            GrB_Index *tr_rows, *tr_cols ;
            uint32_t *tr_vals ;

            OK (LAGraph_Random_Matrix (&A_tr, GrB_UINT32, n, m, tests [k].density, tests [k].seed, msg)) ;

            OK (GrB_Matrix_nvals (&A_tr_nvals, A_tr)) ;

            OK (GrB_Matrix_new (&A, GrB_UINT32, n + m, n + m)) ;

            OK (LAGraph_Malloc ((void**)(&tr_rows), A_tr_nvals, sizeof(GrB_Index), msg)) ;
            OK (LAGraph_Malloc ((void**)(&tr_cols), A_tr_nvals, sizeof(GrB_Index), msg)) ;
            OK (LAGraph_Malloc ((void**)(&tr_vals), A_tr_nvals, sizeof(uint32_t), msg)) ;

            OK (GrB_Matrix_extractTuples (tr_rows, tr_cols, tr_vals, &A_tr_nvals, A_tr)) ;

            for (GrB_Index i = 0; i < A_tr_nvals; i++) {
                GrB_Index row = tr_rows[i];
                GrB_Index col = tr_cols[i];
                uint32_t val = tr_vals[i];
                OK (GrB_Matrix_setElement (A, val, row, col + n)) ;
                OK (GrB_Matrix_setElement (A, val, col + n, row)) ;
            }

            OK (GrB_free (&A_tr)) ;
            OK (LAGraph_Free ((void**)(&tr_rows), msg)) ;
            OK (LAGraph_Free ((void**)(&tr_cols), msg)) ;
            OK (LAGraph_Free ((void**)(&tr_vals), msg)) ;

        } else {
            GrB_Index n = tests [k].n ;

            GrB_Matrix A_dup = NULL ;
            GrB_Index nvals ;
            GrB_Index *rows, *cols ;
            uint32_t *vals ;

            OK (LAGraph_Random_Matrix (&A_dup, GrB_UINT32, n, n, tests [k].density, tests [k].seed, msg)) ;
            OK (GrB_Matrix_new (&A, GrB_UINT32, n, n)) ;

            OK (GrB_Matrix_nvals (&nvals, A_dup)) ;

            OK (LAGraph_Malloc ((void**)(&rows), nvals, sizeof(GrB_Index), msg)) ;
            OK (LAGraph_Malloc ((void**)(&cols), nvals, sizeof(GrB_Index), msg)) ;
            OK (LAGraph_Malloc ((void**)(&vals), nvals, sizeof(uint32_t), msg)) ;

            OK (GrB_Matrix_extractTuples (rows, cols, vals, &nvals, A_dup)) ;

            for (GrB_Index i = 0; i < nvals; i++) {
                GrB_Index row = rows[i];
                GrB_Index col = cols[i];
                uint32_t val = vals[i];
                if (col < row){
                    // use lower triangular entries for the entire matrix
                    OK (GrB_Matrix_setElement (A, val, col, row)) ;
                    OK (GrB_Matrix_setElement (A, val, row, col)) ;
                }
            }
            
            OK (GrB_free (&A_dup)) ;
            OK (LAGraph_Free ((void**)(&rows), msg)) ;
            OK (LAGraph_Free ((void**)(&cols), msg)) ;
            OK (LAGraph_Free ((void**)(&vals), msg)) ;
        }
        // old code using files
        //--------------
        /*
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        TEST_MSG ("Filename %s is invalid", filename) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        */
        //--------------

        TEST_CHECK (A != NULL) ;
        TEST_MSG ("Building of adjacency matrix failed") ;

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

        OK (LAGraph_Incidence_Matrix (&E, G, msg)) ;
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
                // random
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
    OK (LAGraph_Random_Finalize (msg)) ;
}

void test_MaximalMatchingErrors (void)
{
    OK (LAGraph_Init (msg)) ;
//  GrB_set (GrB_GLOBAL, (int32_t) (true), GxB_BURBLE) ;

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

    OK (LAGraph_Finalize (msg)) ;
}

TEST_LIST = {
    { "MaximalMatching", test_MaximalMatching },
    { "MaximalMatchingErrors", test_MaximalMatchingErrors },
    { NULL, NULL }
} ;
