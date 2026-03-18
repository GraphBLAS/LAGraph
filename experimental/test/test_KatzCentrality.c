//------------------------------------------------------------------------------
// LAGraph/src/test/test_KatzCentrality.c: testing for Katz centrality 
// -----------------------------------------------------------------------------

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
//------------------------------------------------------------------------------

double karate_katz [34] = {
    4.9829901679, 3.6518079520, 4.1214049703, 3.0226461989, 1.8904511529,
    2.0310725862, 2.0310725862, 2.5778843937, 3.1126620496, 1.9260736983,
    1.8904511529, 1.4982988513, 1.8005633729, 3.0918177427, 1.9405256847,
    1.9405256847, 1.4062144433, 1.8634795227, 1.9405256847, 2.3774128717,
    1.9405256847, 1.8634795227, 1.9405256847, 2.5865313998, 1.7091149570,
    1.7301052384, 1.7513648061, 2.3556382376, 2.2266143974, 2.3743152002,
    2.6169724552, 3.0054078802, 4.2659247944, 5.1393352269
} ;

double diamonds_katz [8] = {
    1.0000000000, 1.1000000000, 1.3200000000, 1.1000000000, 1.3520000000,
    1.1352000000, 1.1352000000, 1.2270400000
} ;

//------------------------------------------------------------------------------
// difference: compare Katz vector result with reference values
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


void test_katz_diamonds (void)
{
#if LAGRAPH_SUITESPARSE
	LAGraph_Graph G = NULL ;
	OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
	GrB_Vector centrality = NULL ;
    int64_t niters = 0 ;

    // Create diamonds graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "diamonds.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    // Check that AT is cached
    int result = LAGraph_Cached_AT (G, msg) ;
	TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    // Print graph stats
    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows(&n, G->A)) ;
    OK (GrB_Matrix_nvals(&nedges, G->A)) ;
    printf ("\n\nDiamonds graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n", n, nedges) ;

    double alpha = 0.1 ;

    // Compute katz centrality
    double t = LAGraph_WallClockTime() ;
    OK (LAGr_KatzCentrality (&centrality, &niters, G, alpha, 1.0, 1000, 1e-6, false, false, msg)) ;
    t = LAGraph_WallClockTime() - t ;
    printf ("  Time for LAGr_KatzCentrality: %g sec\n", t) ;
    printf ("  Iterations for LAGr_KatzCentrality: %" PRId64 "\n", niters) ;

	// Compare with reference values.
	GrB_Index cn = 0, cnvals = 0 ;
	OK (GrB_Vector_size (&cn, centrality)) ;
	OK (GrB_Vector_nvals (&cnvals, centrality)) ;
	TEST_CHECK (cn == n) ;
	TEST_CHECK (cnvals == n) ;

	double err = difference (centrality, diamonds_katz, 8) ;
	printf ("  diamonds: err: %e\n", err) ;
	TEST_CHECK (err < 1e-4) ;

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
	OK (LAGraph_Finalize (msg)) ;

#endif
}

void test_katz_karate (void)
{
#if LAGRAPH_SUITESPARSE
	LAGraph_Graph G = NULL ;
	OK (LAGraph_Init (msg)) ;
    GrB_Matrix A = NULL ;
	GrB_Vector centrality = NULL ;
    int64_t niters = 0 ;

    // Create karate graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    // Check that AT is cached
    int result = LAGraph_Cached_AT (G, msg) ;
	TEST_CHECK (result == GrB_SUCCESS || result == LAGRAPH_CACHE_NOT_NEEDED) ;

    // Print graph stats
    uint64_t n, nedges ;
    OK (GrB_Matrix_nrows(&n, G->A)) ;
    OK (GrB_Matrix_nvals(&nedges, G->A)) ;
    printf ("\n\nKarate graph (%" PRIu64 " nodes, %" PRIu64 " edges):\n", n, nedges) ;

    // alpha must be less than 1/lambda_max, where lambda_max is 6.725697727631724 for karate graph
    double alpha = 0.1 ;

    // Compute katz centrality
    double t = LAGraph_WallClockTime() ;
    OK (LAGr_KatzCentrality (&centrality, &niters, G, alpha, 1.0, 1000, 1e-6, false, false, msg)) ;
    t = LAGraph_WallClockTime() - t ;
    printf ("  Time for LAGr_KatzCentrality: %g sec\n", t) ;
    printf ("  Iterations for LAGr_KatzCentrality: %" PRId64 "\n", niters) ;

	// Compare with reference values.
	GrB_Index cn = 0, cnvals = 0 ;
	OK (GrB_Vector_size (&cn, centrality)) ;
	OK (GrB_Vector_nvals (&cnvals, centrality)) ;
	TEST_CHECK (cn == n) ;
	TEST_CHECK (cnvals == n) ;

	double err = difference (centrality, karate_katz, 34) ;
	printf ("  karate: err: %e\n", err) ;
	TEST_CHECK (err < 1e-4) ;

    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;
	OK (LAGraph_Finalize (msg)) ;

#endif
}


//------------------------------------------------------------------------------
// list of tests
//------------------------------------------------------------------------------

TEST_LIST = {
    {"test_katz_diamonds", test_katz_diamonds},
	{"test_katz_karate", test_katz_karate},
	{NULL, NULL}
} ;
