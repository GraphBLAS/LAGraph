//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_HarmonicCentrality.c: test cases for
// harmonic centrality (approximate HLL and exact BFS)
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel A. Gomez, FalkorDB

//------------------------------------------------------------------------------

// NOTE: these tests require SuiteSparse:GraphBLAS

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
GrB_Matrix A = NULL ;
GrB_Vector scores_approx = NULL ;
GrB_Vector scores_exact = NULL ;
GrB_Vector node_weights = NULL ;
#define LEN 512
char filename [LEN+1] ;

const char *files [ ] =
{
    "karate.mtx",
    "A.mtx",
    "jagmesh7.mtx",
    "ldbc-directed-example.mtx",
    "ldbc-undirected-example.mtx",
    "",
} ;

//------------------------------------------------------------------------------
// max_relative_error: max percentage difference between two FP64 vectors
//------------------------------------------------------------------------------
// For each entry i where exact[i] != 0, compute |approx[i] - exact[i]| /
// |exact[i]|. Returns the maximum of these ratios.

double max_relative_error
(
    GrB_Vector approx,
    GrB_Vector exact
)
{
    GrB_Index n ;
    GrB_Vector_size (&n, exact) ;

    double max_err = 0 ;
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        double a = 0, e = 0 ;
        GrB_Info info_a = GrB_Vector_extractElement_FP64 (&a, approx, i) ;
        GrB_Info info_e = GrB_Vector_extractElement_FP64 (&e, exact, i) ;

        // skip entries not present in both vectors
        if (info_a != GrB_SUCCESS || info_e != GrB_SUCCESS) continue ;

        if (fabs (e) > 0)
        {
            double rel = fabs (a - e) / fabs (e) ;
            if (rel > max_err) max_err = rel ;
        }
    }
    return max_err ;
}

//------------------------------------------------------------------------------
// test_HarmonicCentrality: compare approximate vs exact on small graphs
//------------------------------------------------------------------------------

void test_HarmonicCentrality (void)
{
    #if LAGRAPH_SUITESPARSE
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; k++)
    {
        const char *aname = files [k] ;
        if (strlen (aname) == 0) break ;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // get matrix dimensions
        GrB_Index n ;
        OK (GrB_Matrix_nrows (&n, A)) ;

        // create boolean node_weights (all nodes, weight = true)
        OK (GrB_Vector_new (&node_weights, GrB_BOOL, n)) ;
        OK (GrB_Vector_assign_BOOL (
            node_weights, NULL, NULL, true, GrB_ALL, n, NULL)) ;

        // compute approximate harmonic centrality
        OK (LAGr_HarmonicCentrality (
            &scores_approx, NULL, A, node_weights, msg)) ;
        printf ("%s", msg);

        // compute exact harmonic centrality
        OK (LAGr_HarmonicCentrality_exact (
            &scores_exact, NULL, A, node_weights, node_weights, msg)) ;

        // print results for small graphs
        LAGraph_PrintLevel pr =
            (n <= 40) ? LAGraph_COMPLETE : LAGraph_SHORT ;
        printf ("\napproximate scores:\n") ;
        OK (LAGraph_Vector_Print (scores_approx, pr, stdout, msg)) ;
        printf ("\nexact scores:\n") ;
        OK (LAGraph_Vector_Print (scores_exact, pr, stdout, msg)) ;

        // compare: max relative (percentage) error
        double rel_err = max_relative_error (scores_approx, scores_exact) ;
        printf ("max relative error: %.2f%%\n", rel_err * 100) ;

        // allow up to 50% relative error for HLL on small graphs
        TEST_CHECK (rel_err < 0.5) ;
        TEST_MSG ("Relative error too large: %.2f%%", rel_err * 100) ;

        // cleanup for this iteration
        OK (GrB_free (&A)) ;
        OK (GrB_free (&scores_approx)) ;
        OK (GrB_free (&scores_exact)) ;
        OK (GrB_free (&node_weights)) ;
    }

    LAGraph_Finalize (msg) ;
    #endif
}

//------------------------------------------------------------------------------
// test_HarmonicCentrality_empty: test with empty node set
//------------------------------------------------------------------------------

void test_HarmonicCentrality_empty (void)
{
    #if LAGRAPH_SUITESPARSE
    LAGraph_Init (msg) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;

    GrB_Index n ;
    OK (GrB_Matrix_nrows (&n, A)) ;

    // empty node_weights — no participating nodes
    OK (GrB_Vector_new (&node_weights, GrB_BOOL, n)) ;

    OK (LAGr_HarmonicCentrality (&scores_approx, NULL, A, node_weights, msg)) ;

    // scores should exist but have no entries
    GrB_Index nvals ;
    OK (GrB_Vector_nvals (&nvals, scores_approx)) ;
    TEST_CHECK (nvals == 0) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&scores_approx)) ;
    OK (GrB_free (&node_weights)) ;

    LAGraph_Finalize (msg) ;
    #endif
}

//------------------------------------------------------------------------------
// test_errors: test error handling
//------------------------------------------------------------------------------

void test_errors (void)
{
    #if LAGRAPH_SUITESPARSE
    LAGraph_Init (msg) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;

    GrB_Index n ;
    OK (GrB_Matrix_nrows (&n, A)) ;
    OK (GrB_Vector_new (&node_weights, GrB_BOOL, n)) ;
    OK (GrB_Vector_assign_BOOL (
        node_weights, NULL, NULL, true, GrB_ALL, n, NULL)) ;

    int result ;

    // scores is NULL
    result = LAGr_HarmonicCentrality (NULL, NULL, A, node_weights, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // A is NULL
    result = LAGr_HarmonicCentrality (
        &scores_approx, NULL, NULL, node_weights, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // node_weights is NULL
    result = LAGr_HarmonicCentrality (
        &scores_approx, NULL, A, NULL, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&node_weights)) ;

    LAGraph_Finalize (msg) ;
    #endif
}

//****************************************************************************

TEST_LIST = {
    {"HarmonicCentrality", test_HarmonicCentrality},
    {"HarmonicCentrality_empty", test_HarmonicCentrality_empty},
    {"HarmonicCentrality_errors", test_errors},
    {NULL, NULL}
} ;
