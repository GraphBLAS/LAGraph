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
GrB_Vector *parent = NULL, *newlabels = NULL ;  // comes from the Coarsen_Matching function
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
    {5, 0.3, 23, LAGraph_Matching_random, 0, 1, "small-random-nopreserve-combine"},
    {0, 0, 0, 0, 0, 0, ""}
} ;

#define LEN 512
#define SEEDS_PER_TEST 1

char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

void test_Coarsen_Matching () {

    OK (LAGraph_Init (msg)) ;
    OK (LAGraph_Random_Init (msg)) ;

    for (int k = 0 ; ; k++)
    {
        const char *aname = tests [k].name ;
        if (strlen (aname) == 0) break ;
        TEST_CASE (aname) ;

        // first generate graph (code from test_MaximalMatching) ===========
        GrB_Index n = tests [k].n ;

        GrB_Matrix A_dup = NULL ;
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
        // graph generation done ======================================

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
                &newlabels,
                G,
                tests [k].matching_type, 
                tests [k].preserve_mapping,
                tests [k].combine_weights,
                1,
                matching_seed,
                msg
            )) ;

            // TODO: Check parent vector for correctness (must be derived from a valid matching)
            OK (LG_check_coarsen (
                &A_coarse_naive,
                G->A,
                parent[0],
                newlabels[0],
                tests [k].preserve_mapping,
                tests [k].combine_weights,
                msg
            )) ;
            OK (LAGraph_Matrix_IsEqual (&ok, A_coarse_LAGraph, A_coarse_naive, msg)) ;
            TEST_CHECK (ok) ;
            printf ("Coarsened matrices do not match for test: %s", tests [k].name) ;
            // OK (LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg)) ;
            // printf("isnull? %d\n", A_coarse_LAGraph == NULL) ;
            OK (LAGraph_Matrix_Print(G->A, LAGraph_COMPLETE, stdout, msg)) ;
            OK (LAGraph_Matrix_Print(A_coarse_LAGraph, LAGraph_COMPLETE, stdout, msg)) ;
            OK (LAGraph_Matrix_Print(A_coarse_naive, LAGraph_COMPLETE, stdout, msg)) ;
            matching_seed += tests [k].n ;
        }
    }
}

TEST_LIST = {
    {"Coarsen_Matching", test_Coarsen_Matching},
    {NULL, NULL}
} ;