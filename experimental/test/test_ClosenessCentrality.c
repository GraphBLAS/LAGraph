//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_ClosenessCentrality.c: tests for Closeness
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Karan Bhalla and Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <acutest.h>

#include "LAGraphX.h"
#include "LAGraph_test.h"
#include "LG_internal.h"

#define LEN 512
char msg [LAGRAPH_MSG_LEN] ;
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// reference results from NetworkX
// note: WF_improved is set to FALSE for all tests 
//------------------------------------------------------------------------------

double diamonds_closeness[8] = {
    0.0, 1.0, 1.0, 1.0, 0.8, 0.5, 0.5, 0.4117647058823529
} ;

double karate_closeness[34] = {
    0.5689655172413793, 0.4852941176470588, 0.559322033898305, 0.4647887323943662,
    0.3793103448275862, 0.38372093023255816, 0.38372093023255816, 0.44,
    0.515625, 0.4342105263157895, 0.3793103448275862, 0.36666666666666664,
    0.3707865168539326, 0.515625, 0.3707865168539326, 0.3707865168539326,
    0.28448275862068967, 0.375, 0.3707865168539326, 0.5,
    0.3707865168539326, 0.375, 0.3707865168539326, 0.39285714285714285,
    0.375, 0.375, 0.3626373626373626, 0.4583333333333333,
    0.4520547945205479, 0.38372093023255816, 0.4583333333333333, 0.5409836065573771,
    0.515625, 0.55
} ;

double ldbc_directed_example_closeness[10] = {
    1.6891891892, 0.0000000000, 1.5151515152, 1.3937282230,
    1.8115942029, 0.0000000000, 0.0000000000, 3.2258064516,
    0.0000000000, 1.1928429423
} ;

//------------------------------------------------------------------------------
// difference: compare closeness vector result with reference values
//------------------------------------------------------------------------------

double difference (GrB_Vector c, double *reference_c, GrB_Index n) ;

double difference (GrB_Vector c, double *reference_c, GrB_Index n)
{
    GrB_Vector diff = NULL, reference_c_vector = NULL ;
    OK (GrB_Vector_new (&reference_c_vector, GrB_FP64, n)) ;

    for (GrB_Index i = 0 ; i < n ; i++)
    {
        OK (GrB_Vector_setElement_FP64 (reference_c_vector, reference_c [i], i)) ;
    }

    OK (GrB_Vector_new (&diff, GrB_FP64, n)) ;
    OK (GrB_eWiseAdd (diff, NULL, NULL, GrB_MINUS_FP64, reference_c_vector, c,
        NULL)) ;
    OK (GrB_apply (diff, NULL, NULL, GrB_ABS_FP64, diff, NULL)) ;

    double err = 0 ;
    OK (GrB_reduce (&err, NULL, GrB_MAX_MONOID_FP64, diff, NULL)) ;

    OK (GrB_free (&diff)) ;
    OK (GrB_free (&reference_c_vector)) ;

    return err ;
}

//------------------------------------------------------------------------------
// directed and unweighted graph, BFS, all sources
//------------------------------------------------------------------------------
void test_closeness_diamonds_bfs (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Graph G = NULL ;
    OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "diamonds.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    OK (GrB_Matrix_nvals (&nedges, G->A)) ;
    printf ("\n\nDiamonds graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n",
            n, nedges) ;

    double t = LAGraph_WallClockTime () ;
    OK (LAGr_ClosenessCentrality (&centrality, G, NULL,
                                  false, CC_BFS, NULL, msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("  Time for LAGr_ClosenessCentrality: %g sec\n", t) ;

    GrB_Index cn = 0, cnvals = 0 ;
    OK (GrB_Vector_size (&cn, centrality)) ;
    OK (GrB_Vector_nvals (&cnvals, centrality)) ;
    TEST_CHECK (cn == n) ;

    double err = difference (centrality, diamonds_closeness, 8) ;
    printf ("  diamonds: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;  

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

//------------------------------------------------------------------------------
// undirected and unweighted graph, BFS, all sources
//------------------------------------------------------------------------------
void test_closeness_karate_bfs (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Graph G = NULL ;
    OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    OK (GrB_Matrix_nvals (&nedges, G->A)) ;
    printf ("\n\nKarate graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n",
            n, nedges) ;

    double t = LAGraph_WallClockTime () ;
    OK (LAGr_ClosenessCentrality (&centrality, G, NULL,
                                  false, CC_BFS, NULL, msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("  Time for LAGr_ClosenessCentrality: %g sec\n", t) ;

    GrB_Index cn = 0, cnvals = 0 ;
    OK (GrB_Vector_size (&cn, centrality)) ;
    OK (GrB_Vector_nvals (&cnvals, centrality)) ;
    TEST_CHECK (cn == n) ;

    double err = difference (centrality, karate_closeness, 34) ;
    printf ("  karate: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;  

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

//------------------------------------------------------------------------------
// undirected and unweighted graph, Floyd-Warshall, all sources
//------------------------------------------------------------------------------
void test_closeness_karate_fw (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Graph G = NULL ;
    OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    OK (GrB_Matrix_nvals (&nedges, G->A)) ;
    printf ("\n\nKarate graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n",
            n, nedges) ;

    double t = LAGraph_WallClockTime () ;
    OK (LAGr_ClosenessCentrality (&centrality, G, NULL,
                                  false, CC_FLOYD_WARSHALL, NULL, msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("  Time for LAGr_ClosenessCentrality: %g sec\n", t) ;
    
    GrB_Index cn = 0, cnvals = 0 ;
    OK (GrB_Vector_size (&cn, centrality)) ;
    OK (GrB_Vector_nvals (&cnvals, centrality)) ;
    TEST_CHECK (cn == n) ;

    double err = difference (centrality, karate_closeness, 34) ;
    printf ("  karate: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;  

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

//------------------------------------------------------------------------------
// directed and weighted graph, Floyd-Warshall, all sources
//------------------------------------------------------------------------------
void test_closeness_ldbc_directed_example_fw (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Graph G = NULL ;
    OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "ldbc-directed-example.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    OK (GrB_Matrix_nvals (&nedges, G->A)) ;
    printf ("\n\nldbc directed example graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n",
            n, nedges) ;

    double t = LAGraph_WallClockTime () ;
    OK (LAGr_ClosenessCentrality (&centrality, G, NULL,
                                  true, CC_FLOYD_WARSHALL, NULL, msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("  Time for LAGr_ClosenessCentrality: %g sec\n", t) ;
    
    GrB_Index cn = 0, cnvals = 0 ;
    OK (GrB_Vector_size (&cn, centrality)) ;
    OK (GrB_Vector_nvals (&cnvals, centrality)) ;
    TEST_CHECK (cn == n) ;

    double err = difference (centrality, ldbc_directed_example_closeness, 10) ;
    printf ("  ldbc directed example: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;  

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

//------------------------------------------------------------------------------
// directed and weighted graph, Bellman-Ford, all sources
//------------------------------------------------------------------------------
void test_closeness_ldbc_directed_example_bf (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Graph G = NULL ;
    OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "ldbc-directed-example.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    OK (GrB_Matrix_nvals (&nedges, G->A)) ;
    printf ("\n\nldbc directed example graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n",
            n, nedges) ;

    double t = LAGraph_WallClockTime () ;
    OK (LAGr_ClosenessCentrality (&centrality, G, NULL,
                                  true, CC_BELLMAN_FORD, NULL, msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("  Time for LAGr_ClosenessCentrality: %g sec\n", t) ;
    
    GrB_Index cn = 0, cnvals = 0 ;
    OK (GrB_Vector_size (&cn, centrality)) ;
    OK (GrB_Vector_nvals (&cnvals, centrality)) ;
    TEST_CHECK (cn == n) ;

    double err = difference (centrality, ldbc_directed_example_closeness, 10) ;
    printf ("  ldbc directed example: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;  

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

//------------------------------------------------------------------------------
// directed and weighted graph, SSSP, all sources
//------------------------------------------------------------------------------
void test_closeness_ldbc_directed_example_sssp (void)
{
#if LAGRAPH_SUITESPARSE
    LAGraph_Graph G = NULL ;
    OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "ldbc-directed-example.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    OK (GrB_Matrix_nvals (&nedges, G->A)) ;
    printf ("\n\nldbc directed example graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n",
            n, nedges) ;

    GrB_Scalar Delta = NULL ;
    OK (GrB_Scalar_new (&Delta, GrB_FP64)) ;
    OK (GrB_Scalar_setElement_FP64 (Delta, 1.6)) ;

    double t = LAGraph_WallClockTime () ;
    OK (LAGr_ClosenessCentrality (&centrality, G, NULL,
                                  true, CC_SSSP, Delta, msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("  Time for LAGr_ClosenessCentrality: %g sec\n", t) ;
    
    GrB_Index cn = 0, cnvals = 0 ;
    OK (GrB_Vector_size (&cn, centrality)) ;
    OK (GrB_Vector_nvals (&cnvals, centrality)) ;
    TEST_CHECK (cn == n) ;

    double err = difference (centrality, ldbc_directed_example_closeness, 10) ;
    printf ("  ldbc directed example: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;  

    OK (GrB_free (&Delta)) ;
    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
    OK (LAGraph_Finalize (msg)) ;
#endif
}

//------------------------------------------------------------------------------
// list of tests    
//------------------------------------------------------------------------------

TEST_LIST = {
    {"test_closeness_diamonds_bfs", test_closeness_diamonds_bfs},
    {"test_closeness_karate_bfs",   test_closeness_karate_bfs},
    {"test_closeness_karate_fw",    test_closeness_karate_fw},
    {"test_closeness_ldbc_directed_example_fw", test_closeness_ldbc_directed_example_fw},
    {"test_closeness_ldbc_directed_example_bf", test_closeness_ldbc_directed_example_bf},
    {"test_closeness_ldbc_directed_example_sssp", test_closeness_ldbc_directed_example_sssp},
    {NULL, NULL}
} ;
