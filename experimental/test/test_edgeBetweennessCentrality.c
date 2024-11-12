//------------------------------------------------------------------------------
// LAGraph/src/test/test_Betweenness.c: test cases for BC (GAP method)
// -----------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>

#include <LAGraph_test.h>

#define LEN 512
char msg [LAGRAPH_MSG_LEN] ;
char filename [LEN+1] ;
LAGraph_Graph G = NULL ;

//------------------------------------------------------------------------------
// difference: compare the LAGraph and GAP results
//------------------------------------------------------------------------------

float difference (GrB_Vector bc, double *gap_result) ;

float difference (GrB_Vector bc, double *gap_result)
{
    GrB_Vector diff = NULL, gap_bc = NULL ;
    GrB_Index n = 0 ;
    OK (GrB_Vector_size (&n, bc)) ;
    OK (GrB_Vector_new (&gap_bc, GrB_FP32, n)) ;
    for (int i = 0 ; i < n ; i++)
    {
        OK (GrB_Vector_setElement_FP64 (gap_bc, gap_result [i], i)) ;
    }
    // diff = max (abs (gap_bc - bc))
    OK (GrB_Vector_new (&diff, GrB_FP32, n)) ;
    OK (GrB_eWiseAdd (diff, NULL, NULL, GrB_MINUS_FP32, gap_bc, bc,
        NULL)) ;
    OK (GrB_apply (diff, NULL, NULL, GrB_ABS_FP32, diff, NULL)) ;
    float err = 0 ;
    OK (GrB_reduce (&err, NULL, GrB_MAX_MONOID_FP32, diff, NULL)) ;
    OK (GrB_free (&diff)) ;
    OK (GrB_free (&gap_bc)) ;
    return (err) ;
}


//------------------------------------------------------------------------------
// results for book graph
//------------------------------------------------------------------------------

// Exact results from edge_betweenness_centrality from NetworkX of the graph
// from the book

double example_bc [11] = {
    0.0952380952380952
    0.08333333333333331
    0.0952380952380952
    0.047619047619047616
    0.20238095238095238
    0.19047619047619047
    0.047619047619047616
    0.20238095238095238
    0.2857142857142857
    0.2857142857142857
    0.14285714285714285
    0.14285714285714285 } ; 

//------------------------------------------------------------------------------
// results for karate graph
//------------------------------------------------------------------------------

// Exact results from edge_betweenness_centrality from NetworkX of the karate 
// graph

double karate_bc [78] = {
    0.025252525252525245,
    0.0777876807288572,
    0.02049910873440285,
    0.0522875816993464,
    0.07813428401663694,
    0.07813428401663695,
    0.0228206434088787,
    0.07423959482783014,
    0.0522875816993464,
    0.058823529411764705,
    0.04652406417112298,
    0.04237189825425121,
    0.04012392835922248,
    0.045936960642843,
    0.040123928359222474,
    0.1272599949070537,
    0.023232323232323233,
    0.0077243018419489,
    0.007422969187675069,
    0.01240556828792123,
    0.01869960105254222,
    0.014633732280791102,
    0.01869960105254222,
    0.032280791104320514,
    0.022430184194890075,
    0.025214328155504617,
    0.009175791528732704,
    0.030803836686189627,
    0.007630931160342923,
    0.04119203236850296,
    0.02278244631185807,
    0.06898678663384543,
    0.003365588659706307,
    0.012299465240641705,
    0.01492233256939139,
    0.0047534165181224,
    0.0029708853238265,
    0.0029708853238265003,
    0.0047534165181224,
    0.029411764705882353,
    0.029411764705882353,
    0.00980392156862745,
    0.0304416716181422,
    0.04043657867187279,
    0.029615482556659026,
    0.06782389723566191,
    0.024083977025153497,
    0.03473955238661121,
    0.024083977025153497,
    0.03473955238661121,
    0.024083977025153497,
    0.03473955238661121,
    0.05938233879410351,
    0.024083977025153497,
    0.03473955238661121,
    0.024083977025153493,
    0.03473955238661121,
    0.019776193305605066,
    0.010536739948504653,
    0.00665478312537136,
    0.022341057635175278,
    0.03266983561101209,
    0.0042186571598336305,
    0.018657159833630418,
    0.040106951871657755,
    0.04205783323430383,
    0.004532722179781003,
    0.0542908072319837,
    0.030477039300568713,
    0.0148544266191325,
    0.024564977506153975,
    0.023328523328523323,
    0.029807882749059215,
    0.01705288175876411,
    0.02681436210847975,
    0.04143394731630026,
    0.05339388280564752,
    0.008225108225108224 } ; 


//------------------------------------------------------------------------------
// test_bc
//------------------------------------------------------------------------------

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
    OK (LG_check_edgeBetweennessCentrality (&centrality, G, msg)) ;
    printf ("\nkarate bc:\n") ;
    OK (LAGraph_Delete (&G, msg)) ;

    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_bc_brutal: test BetweenessCentraliy with brutal malloc debugging
//------------------------------------------------------------------------------

#if LAGRAPH_SUITESPARSE
void test_bc_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

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
    printf ("\n") ;

    // compute its betweenness centrality
    LG_BRUTAL_BURBLE (LAGr_Betweenness (&centrality, G,
            karate_sources, 4, msg)) ;

    // compare with GAP:
    float err = difference (centrality, karate_bc) ;
    printf ("karate:   err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//------------------------------------------------------------------------------
// list of tests
//------------------------------------------------------------------------------

TEST_LIST = {
    {"test_bc", test_bc},
    #if LAGRAPH_SUITESPARSE
    {"test_bc_brutal", test_bc_brutal },
    #endif
    {NULL, NULL}
};
