//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_CoarsenMatching
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
#include "LAGraphX.h"
#include "LAGraph_test.h"
#include "LG_Xtest.h"
#include "LG_internal.h"

char msg [LAGRAPH_MSG_LEN] ;

GrB_Matrix A = NULL, A_coarse_LAGraph = NULL, A_coarse_naive = NULL ;
GrB_Vector parent = NULL, newlabel = NULL, inv_newlabel = NULL ;    // outputs from the Coarsen_Matching function
LAGraph_Graph G = NULL ;

typedef struct
{
    // options for building graph:
    GrB_Index n ;         // number of nodes in the graph
    double density ;      // density of the matrix
    uint64_t seed ;       // seed used to generate the graph for this test
    
    // options for coarsening (see Coarsen_Matching function for details):
    LAGraph_Matching_kind matching_type ;
    int preserve_mapping ;
    int combine_weights ;

    const char *name ;
}
matrix_info ;

const matrix_info tests [ ] = {
    // random, preserve, combine
    {10, 0.3, 55, LAGraph_Matching_random, 1, 1, "small-random-preserve-combine"},
    {500, 0.4, 16, LAGraph_Matching_random, 1, 1, "large-random-preserve-combine"},
    // random, preserve, nocombine
    {10, 0.3, 62, LAGraph_Matching_random, 1, 0, "small-random-preserve-nocombine"},
    {500, 0.4, 21, LAGraph_Matching_random, 1, 0, "large-random-preserve-nocombine"},
    // random, nopreserve, combine
    {10, 0.3, 23, LAGraph_Matching_random, 0, 1, "small-random-nopreserve-combine"},
    {500, 0.4, 31, LAGraph_Matching_random, 0, 1, "large-random-nopreserve-combine"},
    // random, nopreserve, nocombine
    {10, 0.3, 92, LAGraph_Matching_random, 0, 0, "small-random-nopreserve-nocombine"},
    {500, 0.4, 44, LAGraph_Matching_random, 0, 0, "large-random-nopreserve-nocombine"},

    // same as above except weighted matching (mix of light and heavy)
    // random, preserve, combine
    {10, 0.3, 55, LAGraph_Matching_heavy, 1, 1, "small-random-preserve-combine"},
    {500, 0.4, 16, LAGraph_Matching_light, 1, 1, "large-random-preserve-combine"},
    // random, preserve, nocombine
    {10, 0.3, 62, LAGraph_Matching_light, 1, 0, "small-random-preserve-nocombine"},
    {500, 0.4, 21, LAGraph_Matching_heavy, 1, 0, "large-random-preserve-nocombine"},
    // random, nopreserve, combine
    {10, 0.3, 23, LAGraph_Matching_light, 0, 1, "small-random-nopreserve-combine"},
    {500, 0.4, 31, LAGraph_Matching_heavy, 0, 1, "large-random-nopreserve-combine"},
    // random, nopreserve, nocombine
    {10, 0.3, 92, LAGraph_Matching_heavy, 0, 0, "small-random-nopreserve-nocombine"},
    {500, 0.4, 44, LAGraph_Matching_light, 0, 0, "large-random-nopreserve-nocombine"},

    {0, 0, 0, 0, 0, 0, ""}
} ;

#define LEN 512
#define SEEDS_PER_TEST 3

char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

void test_Coarsen_Matching () {

    OK (LAGraph_Init (msg)) ;
    OK (LAGraph_Random_Init (msg)) ;

#if LAGRAPH_SUITESPARSE
    for (int k = 0 ; ; k++)
    {
        const char *aname = tests [k].name ;
        if (strlen (aname) == 0) { break ; }
        TEST_CASE (aname) ;

        // ================= first generate graph (code from test_MaximalMatching) =================
        GrB_Index n = tests [k].n ;

        GrB_Matrix A_dup = NULL ;
        A = NULL ;
        GrB_Index nvals ;
        GrB_Index *rows, *cols ;
        double *vals ;

        OK (LAGraph_Random_Matrix (&A_dup, GrB_FP64, n, n, tests [k].density, tests [k].seed, msg)) ;
        OK (GrB_Matrix_new (&A, GrB_FP64, n, n)) ;

        OK (GrB_Matrix_nvals (&nvals, A_dup)) ;

        OK (LAGraph_Malloc ((void**)(&rows), nvals, sizeof(GrB_Index), msg)) ;
        OK (LAGraph_Malloc ((void**)(&cols), nvals, sizeof(GrB_Index), msg)) ;
        OK (LAGraph_Malloc ((void**)(&vals), nvals, sizeof(double), msg)) ;

        OK (GrB_Matrix_extractTuples (rows, cols, vals, &nvals, A_dup)) ;

        for (GrB_Index i = 0; i < nvals; i++) {
            GrB_Index row = rows [i];
            GrB_Index col = cols [i];
            double val = vals [i];
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
        // =============================== graph generation done ======================================

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

        bool ok = 0;
        OK (LAGraph_Matrix_IsEqual (&ok, G->A, G->AT, msg)) ;
        TEST_CHECK (ok) ;
        TEST_MSG ("Input graph is not undirected") ;
        G->kind = LAGraph_ADJACENCY_UNDIRECTED ;

        uint64_t matching_seed = 0 ;
        for (int i = 0; i < SEEDS_PER_TEST ; i++) {
            OK (LAGraph_Coarsen_Matching (
                &A_coarse_LAGraph,
                &parent,
                &newlabel,
                &inv_newlabel,
                G,
                tests [k].matching_type, 
                tests [k].preserve_mapping,
                tests [k].combine_weights,
                matching_seed,
                msg
            )) ;
            OK (LG_check_coarsen (
                &A_coarse_naive,
                G->A,
                parent,
                (newlabel == NULL ? NULL : newlabel),
                (inv_newlabel == NULL ? NULL : inv_newlabel),
                tests [k].preserve_mapping,
                tests [k].combine_weights,
                msg
            )) ;

            if (newlabel != NULL) {
                GrB_free (&newlabel) ;
                GrB_free (&inv_newlabel) ;
            }
            // Check parent vector for matching-specific correctness (must be derived from a valid matching)
            // requirements: no node is the parent of more than 2 nodes, and if p[i] != i, then A[i][p[i]] exists
            int8_t *freq ;
            OK (LAGraph_Malloc ((void**)(&freq), n, sizeof(int8_t), msg)) ;
            memset(freq, 0, n * sizeof(int8_t)) ;

            for (GrB_Index i = 0 ; i < n ; i++) {
                uint64_t par ;
                OK (GrB_Vector_extractElement (&par, parent, i)) ;
                freq [par]++ ;
                TEST_CHECK (freq [par] <= 2) ;
                TEST_MSG ("Parent vector not from a valid matching for test: %s\n", tests [k].name) ;

                if (par != i) {
                    // make sure that (i, par) is a edge in the graph
                    TEST_CHECK (GxB_Matrix_isStoredElement (G->A, i, par) == GrB_SUCCESS) ;
                    TEST_MSG ("Parent vector not from a valid matching for test: %s\n", tests [k].name) ;
                }
            }
            GrB_free (&parent) ;
            
            OK (LAGraph_Free ((void**)(&freq), msg)) ;

#if 0
            OK (LAGraph_Matrix_Print (G->A, LAGraph_COMPLETE, stdout, msg)) ;
            OK (LAGraph_Matrix_Print (A_coarse_LAGraph, LAGraph_COMPLETE, stdout, msg)) ;
            OK (LAGraph_Matrix_Print (A_coarse_naive, LAGraph_COMPLETE, stdout, msg)) ;
            OK (LAGraph_Vector_Print (parent[0], LAGraph_COMPLETE, stdout, msg)) ;
            OK (LAGraph_Vector_Print (newlabels[0], LAGraph_COMPLETE, stdout, msg)) ;
#endif
            OK (LAGraph_Matrix_IsEqual (&ok, A_coarse_LAGraph, A_coarse_naive, msg)) ;
            TEST_CHECK (ok) ;
            TEST_MSG ("Coarsened matrices do not match for test: %s", tests [k].name) ;

            OK (GrB_free (&A_coarse_LAGraph)) ;
            OK (GrB_free (&A_coarse_naive)) ;

            matching_seed += tests [k].n ;
        }
        OK (LAGraph_Delete (&G, msg)) ;
    }
#endif

    OK (LAGraph_Finalize (msg)) ;
    OK (LAGraph_Random_Finalize (msg)) ;
}

void test_Coarsen_Matching_Errors() {

    OK (LAGraph_Init (msg)) ;

#if LAGRAPH_SUITESPARSE
    OK (GrB_Matrix_new (&A, GrB_FP64, 5, 5)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;

    G->kind = LAGraph_ADJACENCY_DIRECTED ;

    // directed graph
    GrB_Info result = LAGraph_Coarsen_Matching (NULL, NULL, NULL, NULL, G, 0, 0, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == LAGRAPH_INVALID_GRAPH) ;

    G->kind = LAGraph_ADJACENCY_UNDIRECTED ;
    G->nself_edges = 1 ;

    // non-zero self-loops
    result = LAGraph_Coarsen_Matching (NULL, NULL, NULL, NULL, G, 0, 0, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == LAGRAPH_NO_SELF_EDGES_ALLOWED) ;

    G->nself_edges = 0;

    // output coarsened matrix pointer is NULL
    result = LAGraph_Coarsen_Matching (NULL, NULL, NULL, NULL, G, 0, 0, 0, 0, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

#endif

    OK (LAGraph_Finalize (msg)) ;
}


void test_Coarsen_Matching_NullInputs() {

    OK (LAGraph_Init (msg)) ;
    OK (LAGraph_Random_Init (msg)) ;

#if LAGRAPH_SUITESPARSE

    OK (GrB_Matrix_new (&A, GrB_FP64, 5, 5)) ;

    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    OK (LAGraph_Cached_NSelfEdges (G, msg)) ;
    OK (LAGraph_Cached_AT (G, msg) < 0) ; // warning is expected; check for error

    // do this to get full code coverage and catch any unexpected behavior
    OK (LAGraph_Coarsen_Matching (&A_coarse_LAGraph, NULL, NULL, NULL, G, 0, 0, 1, 42, msg)) ;

    OK (GrB_free (&A_coarse_LAGraph)) ;
    // do it with parent, inv_newlabels NULL but newlabels not NULL
    OK (LAGraph_Coarsen_Matching (&A_coarse_LAGraph, NULL, &newlabel, NULL, G, 0, 0, 1, 42, msg)) ;

    OK (GrB_free (&A_coarse_LAGraph)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    TEST_CHECK (newlabel != NULL) ;
    TEST_MSG ("Null input check failed!\n") ;

    OK (GrB_free (&newlabel)) ;
#endif

    OK (LAGraph_Finalize (msg)) ;
    OK (LAGraph_Random_Finalize (msg)) ;
}

TEST_LIST = {
    {"Coarsen_Matching", test_Coarsen_Matching},
    {"Coarsen_Matching_Errors", test_Coarsen_Matching_Errors},
    {"Coarsen_Matching_NullInputs", test_Coarsen_Matching_NullInputs},
    {NULL, NULL}
} ;