//------------------------------------------------------------------------------
// LAGraph/src/test/test_PageRankGAP.c: test cases for pagerank (GAP method)
// -----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// brutal HERE TODO FIXME

#include <stdio.h>
#include <acutest.h>

#include <LAGraph_test.h>

#define LEN 512
char msg [LAGRAPH_MSG_LEN] ;
char filename [LEN+1] ;
LAGraph_Graph G = NULL ;

//------------------------------------------------------------------------------
// difference: compare the LAGraph and MATLAB results
//------------------------------------------------------------------------------

float difference (GrB_Vector centrality, double *matlab_result) ;

float difference (GrB_Vector centrality, double *matlab_result)
{
    GrB_Vector diff = NULL, cmatlab = NULL ;
    GrB_Index n = 0 ;
    OK (GrB_Vector_size (&n, centrality)) ;
    OK (GrB_Vector_new (&cmatlab, GrB_FP32, n)) ;
    for (int i = 0 ; i < n ; i++)
    {
        OK (GrB_Vector_setElement_FP64 (cmatlab, matlab_result [i], i)) ;
    }
    // diff = max (abs (cmatlab - centrality))
    OK (GrB_Vector_new (&diff, GrB_FP32, n)) ;
    OK (GrB_eWiseAdd (diff, NULL, NULL, GrB_MINUS_FP32, cmatlab, centrality,
        NULL)) ;
    OK (GrB_apply (diff, NULL, NULL, GrB_ABS_FP32, diff, NULL)) ;
    float err = 0 ;
    OK (GrB_reduce (&err, NULL, GrB_MAX_MONOID_FP32, diff, NULL)) ;
    OK (GrB_free (&diff)) ;
    OK (GrB_free (&cmatlab)) ;
    return (err) ;
}

//------------------------------------------------------------------------------
// valid results for karate graph and west0067 graphs
//------------------------------------------------------------------------------

// These two matrices have no sinks (nodes with zero outdegree) so the MATLAB
// centrality (G, 'pagerank') and LAGraph PageRankGAP results will be the same.

double karate_rank [34] = {
    0.0970011147,
    0.0528720584,
    0.0570750515,
    0.0358615175,
    0.0219857202,
    0.0291233505,
    0.0291233505,
    0.0244945048,
    0.0297681451,
    0.0143104668,
    0.0219857202,
    0.0095668739,
    0.0146475355,
    0.0295415677,
    0.0145381625,
    0.0145381625,
    0.0167900065,
    0.0145622041,
    0.0145381625,
    0.0196092670,
    0.0145381625,
    0.0145622041,
    0.0145381625,
    0.0315206825,
    0.0210719482,
    0.0210013837,
    0.0150430281,
    0.0256382216,
    0.0195723309,
    0.0262863139,
    0.0245921424,
    0.0371606178,
    0.0716632142,
    0.1008786453 } ;

double west0067_rank [67] = {
    0.0233753869,
    0.0139102552,
    0.0123441027,
    0.0145657095,
    0.0142018541,
    0.0100791606,
    0.0128753395,
    0.0143945684,
    0.0110203141,
    0.0110525383,
    0.0119311961,
    0.0072382247,
    0.0188680398,
    0.0141596605,
    0.0174877889,
    0.0170362099,
    0.0120433909,
    0.0219844489,
    0.0195274443,
    0.0394465722,
    0.0112038726,
    0.0090174094,
    0.0140088120,
    0.0122532937,
    0.0153346283,
    0.0135241334,
    0.0158714693,
    0.0149689529,
    0.0144097230,
    0.0137583019,
    0.0314386080,
    0.0092857745,
    0.0081814168,
    0.0102137827,
    0.0096547214,
    0.0129622400,
    0.0244173417,
    0.0173963657,
    0.0127705717,
    0.0143297446,
    0.0140509341,
    0.0104117131,
    0.0173516407,
    0.0149175105,
    0.0119979624,
    0.0095043613,
    0.0153295328,
    0.0077710930,
    0.0259969472,
    0.0126926269,
    0.0088870166,
    0.0080836101,
    0.0096023576,
    0.0091000837,
    0.0246131958,
    0.0159589365,
    0.0183500031,
    0.0155811507,
    0.0157693756,
    0.0116319823,
    0.0230649292,
    0.0149070613,
    0.0157469640,
    0.0134396036,
    0.0189218603,
    0.0114528518,
    0.0223213267 } ;

//------------------------------------------------------------------------------
// tesk_ranker
//------------------------------------------------------------------------------

void test_ranker(void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Type atype = NULL ;
    GrB_Vector centrality = NULL, cmatlab = NULL, diff = NULL ;
    int niters = 0 ;

    // create the karate graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A
    OK (LAGraph_Property_RowDegree (G, msg)) ;

    // compute its pagerank
    OK (LAGraph_VertexCentrality_PageRankGAP (&centrality, G, 0.85,
        1e-4, 100, &niters, msg)) ;
    // OK (LAGraph_Vector_print (centrality, 5, stdout, msg)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    // compare with MATLAB: cmatlab = centrality (G, 'pagerank')
    float err = difference (centrality, karate_rank) ;
    printf ("\nkarate:   err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    // create the west0067 graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "west0067.mtx") ;
    f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A
    OK (LAGraph_Property_AT (G, msg)) ;
    OK (LAGraph_Property_RowDegree (G, msg)) ;

    // compute its pagerank
    OK (LAGraph_VertexCentrality_PageRankGAP (&centrality, G, 0.85,
        1e-4, 100, &niters, msg)) ;
    // OK (LAGraph_Vector_print (centrality, 5, stdout, msg)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    // compare with MATLAB: cmatlab = centrality (G, 'pagerank')
    err = difference (centrality, west0067_rank) ;
    printf ("west0067: err: %e\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// list of tests
//------------------------------------------------------------------------------

TEST_LIST = {
    {"test_ranker", test_ranker},
    {NULL, NULL}
};

