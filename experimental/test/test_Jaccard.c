#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <graph_zachary_karate.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix JC = NULL ;

#define LEN 512
char filename [LEN+1] ;

float A_jaccard_allpairs [21] = 
{
	0.142857, 0.333333, 0.6,      0.5,      0.75,     0.142857,  // (0,1..6)
    0.428571, 0.428571, 0.142857, 0.285714, 1.0,                 // (1,2..6)
    0.428571, 0.6,      0.5,      0.428571,                      // (2,3..6)
    0.333333, 0.5,      0.428571,                                // (3,4..6)
    0.75,     0.142857,                                           // (4,5..6)
    0.285714                                                        // (5,6)
} ;

float karate_jaccard_allpairs[561] = {
    0.388889f, 0.238095f, 0.294118f, 0.117647f, 0.111111f, 0.111111f, 0.176471f,
    0.050000f, 0.117647f, 0.000000f, 0.058824f, 0.166667f, 0.058824f, 0.055556f,
    0.058824f, 0.000000f, 0.000000f, 0.125000f, 0.000000f, 0.000000f, 0.000000f,
    0.000000f, 0.000000f, 0.055556f, 0.055556f, 0.000000f, 0.052632f, 0.117647f,
    0.000000f, 0.111111f, 0.000000f, 0.000000f, 0.000000f, 0.090909f, 0.083333f,
    0.083333f, 0.272727f, 0.100000f, 0.090909f, 0.111111f, 0.222222f, 0.000000f,
    0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
    0.000000f, 0.000000f, 0.000000f, 0.083333f, 0.090909f, 0.000000f, 0.071429f,
    0.105263f, 0.130435f, 0.266667f, 0.363636f, 0.300000f, 0.272727f, 0.100000f,
    0.090909f, 0.100000f, 0.000000f, 0.333333f, 0.272727f, 0.153846f, 0.000000f,
    0.250000f, 0.000000f, 0.000000f, 0.047619f, 0.125000f, 0.111111f, 0.111111f,
    0.222222f, 0.142857f, 0.125000f, 0.166667f, 0.000000f, 0.000000f, 0.333333f,
    0.000000f, 0.285714f, 0.000000f, 0.333333f, 0.000000f, 0.000000f, 0.000000f,
    0.000000f, 0.000000f, 0.000000f, 0.111111f, 0.125000f, 0.000000f, 0.058824f,
    0.045455f, 0.750000f, 0.166667f, 0.142857f, 0.000000f, 0.333333f, 0.250000f,
    0.142857f, 0.000000f, 0.000000f, 0.250000f, 0.250000f, 0.000000f, 0.200000f,
    0.000000f, 0.250000f, 0.000000f, 0.000000f, 0.125000f, 0.000000f, 0.000000f,
    0.142857f, 0.125000f, 0.000000f, 0.250000f, 0.000000f, 0.000000f, 0.000000f,
    0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.200000f, 0.250000f, 0.142857f,
    0.142857f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.200000f, 0.000000f,
    0.200000f, 0.000000f, 0.333333f, 0.200000f, 0.200000f, 0.000000f, 0.000000f,
    0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.111111f, 0.000000f, 0.000000f,
    0.142857f, 0.000000f, 0.000000f, 0.166667f, 0.250000f, 0.125000f, 0.125000f,
    0.166667f, 0.200000f, 0.000000f, 0.200000f, 0.000000f, 0.333333f, 0.166667f,
    0.166667f, 0.000000f, 0.000000f, 0.111111f, 0.200000f, 0.125000f, 0.200000f,
    0.111111f, 0.250000f, 0.088000f, 0.375000f, 0.400000f, 0.142857f, 0.166667f,
    0.100000f, 0.428571f, 0.333333f, 0.400000f, 0.166667f, 0.333333f, 0.400000f,
    0.250000f, 0.000000f, 0.000000f, 0.166667f, 0.285714f, 0.333333f, 0.285714f,
    0.375000f, 0.000000f, 0.000000f, 0.333333f, 0.250000f, 0.333333f, 0.166667f,
    0.200000f, 0.285714f, 0.333333f, 0.200000f, 0.285714f, 0.222222f, 0.200000f,
    0.250000f, 0.400000f, 0.166667f, 0.166667f, 0.111111f, 0.250000f, 0.333333f,
    0.000000f, 0.000000f, 0.400000f, 0.333333f, 0.250000f, 0.000000f, 0.400000f,
    0.000000f, 0.400000f, 0.142857f, 0.200000f, 0.166667f, 0.333333f, 0.400000f,
    0.000000f, 0.000000f, 0.400000f, 0.250000f, 0.400000f, 0.166667f, 0.250000f,
    0.000000f, 0.200000f, 0.000000f, 0.000000f, 0.125000f, 0.142857f, 0.111111f,
    0.111111f, 0.142857f, 0.100000f, 0.000000f, 0.200000f, 0.000000f, 0.000000f,
    0.111111f, 0.000000f, 0.000000f, 0.071429f, 0.000000f, 0.000000f, 0.000000f,
    0.000000f, 0.200000f, 0.250000f, 0.200000f, 0.200000f, 0.000000f, 0.200000f,
    0.000000f, 0.333333f, 0.200000f, 0.285714f, 0.000000f, 0.200000f, 0.125000f,
    0.250000f, 0.000000f, 0.333333f, 0.142857f, 0.090909f, 0.166667f, 0.133333f,
    0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.666667f, 0.000000f, 0.000000f,
    0.142857f, 0.000000f, 0.000000f, 0.333333f, 0.000000f, 0.200000f, 0.285714f,
    0.250000f, 0.200000f, 0.333333f, 0.000000f, 0.200000f, 0.333333f, 0.400000f,
    0.166667f, 0.200000f, 0.166667f, 0.250000f, 0.200000f, 0.285714f, 0.285714f,
    0.333333f, 0.166667f, 0.285714f, 0.333333f, 0.166667f, 0.333333f, 0.285714f,
    0.400000f, 0.166667f, 0.285714f, 0.333333f, 0.125000f, 0.250000f, 0.125000f,
    0.125000f, 0.200000f, 0.200000f, 0.000000f, 0.000000f, 0.333333f, 0.333333f,
    0.333333f, 0.142857f, 0.200000f, 0.166667f, 0.333333f, 0.200000f, 0.200000f,
    0.111111f, 0.111111f, 0.125000f, 0.166667f, 0.133333f, 0.000000f, 0.333333f,
    0.111111f, 0.000000f, 0.111111f, 0.142857f, 0.095238f, 0.526316f,
};

float cover_jaccard_allpairs[21] = {
    0.0f, 0.250000f, 0.0f, 0.250000f, 0.0f, 0.500000f,
    0.2f, 0.5f, 0.200000f, 0.25f, 0.166667f,
    0.200000f, 0.500000f, 0.500000f, 0.166667, 
	0.200000f, 0.250000f, 0.166667f, 
	0.0f, 0.166667f, 
    0.500000f
};



float matrix_difference(GrB_Matrix JC, const float *ref, GrB_Matrix A_mask, bool all_pairs) 
{
    GrB_Matrix JCref = NULL, diff = NULL ;
	GrB_Index N = 0;
	OK (GrB_Matrix_nrows (&N, JC)) ;
    OK (GrB_Matrix_new (&JCref, GrB_FP32, N, N)) ;
	GrB_Index k = 0 ;
    for (GrB_Index i = 0 ; i < N ; i++) {
        for (GrB_Index j = i+1 ; j < N ; j++) {
            float v = ref[k++] ;
            if (v != 0.0f) {
                OK (GrB_Matrix_setElement_FP32 (JCref, v, i, j)) ;
            }
        }
    }
    OK (GrB_Matrix_new (&diff, GrB_FP32, N, N)) ;
	if (all_pairs)
    {
        OK (GrB_eWiseAdd ( diff, NULL, NULL, GrB_MINUS_FP32, JC, JCref, NULL)) ;
    }
    else
    {
        OK (GrB_eWiseAdd ( diff, A_mask, NULL, GrB_MINUS_FP32, JC, JCref, NULL)) ;
    }
    OK (GrB_apply (diff, NULL, NULL, GrB_ABS_FP32, diff, NULL)) ; //absolute value
    float err = 0 ;
    OK (GrB_reduce (&err, NULL, GrB_MAX_MONOID_FP32, diff, NULL)) ; //get max
    OK (GrB_free (&diff)) ;
    OK (GrB_free (&JCref)) ;
    return err ;
}


void run_jaccard (const char *input_mat, const float *results)
{

    GrB_Matrix A = NULL, ABool = NULL;
	GrB_Index n; 

    // create the graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", input_mat) ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;

	// Convert matrix to boolean, symmetrize, and remove self-loops then build graph
	OK( GrB_Matrix_nrows(&n, A));
	OK( GrB_Matrix_new(&ABool, GrB_BOOL, n, n));
	OK( GrB_apply(ABool, NULL, NULL, GrB_ONEB_BOOL, A,(bool) true, NULL));  
	OK( GrB_eWiseAdd(ABool, NULL, NULL, GrB_LOR, ABool, ABool, GrB_DESC_T1));
    OK( GrB_select (ABool, NULL, NULL, GrB_OFFDIAG, ABool, 0, NULL)) ;
    OK( LAGraph_New (&G, &ABool, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (ABool == NULL) ;    

	//--------------------------------------------------------------------------
    // Experiment 1: All-pairs Jaccard similarity
    //--------------------------------------------------------------------------
	printf("\nExp 1 : Jaccard similarity\n");
	bool all_pairs = true;	// if false it only computes jaccard weights(jaccard similarity for neighbors)
    OK (LAGraph_Jaccard(&JC, G, all_pairs, msg)) ;
    FILE *fo = fopen ("LAGraph_Jaccard_all_pairs.txt", "w") ;
	OK (GxB_Matrix_fprint(JC, "my matrix", LAGraph_COMPLETE, fo));
	fclose(fo);

	//check for correctness
	GrB_Matrix A_mask = G->A;
	float err1 = matrix_difference (JC, results, A_mask, all_pairs) ;
    printf("Exp 1 : all-pairs Jaccard, max abs error = %en",err1) ;
    TEST_CHECK (err1 < 1e-4) ;   // adjust tolerance if needed
	OK(GrB_free(&JC)); 
	JC = NULL;

	//--------------------------------------------------------------------------
    // Experiment 2: Jaccard weights (neighbors only)
    //--------------------------------------------------------------------------
	printf("Exp 2 : Jaccard weights\n");
	all_pairs = false;	// if false it only computes jaccard weights(jaccard similarity for neighbors)
    OK (LAGraph_Jaccard(&JC, G, all_pairs, msg)) ;
	// fo = fopen ("LAGraph_Jaccard_weights.txt", "w") ;
	// OK (GxB_Matrix_fprint(JC, "my matrix", LAGraph_COMPLETE, fo));
	// fclose(fo);
   	//check for correctness
	err1 = matrix_difference (JC, results, A_mask, all_pairs) ;
    printf("Exp 1 : all-pairs Jaccard, max abs error = %en",err1) ;
    TEST_CHECK (err1 < 1e-4) ;   // adjust tolerance if needed
    TEST_CHECK(0 == GrB_free(&JC));
	
    // free everything
	OK(LAGraph_Delete(&G, msg));
    OK(GrB_free(&A));
}

void test_Jaccard (void)
{

    LAGraph_Init (msg) ;

	run_jaccard ("A.mtx", A_jaccard_allpairs) ;

	// run_jaccard ("karate.mtx", karate_jaccard_allpairs) ;

	run_jaccard ("cover.mtx", cover_jaccard_allpairs) ;

    OK(LAGraph_Finalize(msg));

}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"test my jaccard", test_Jaccard},    // just one test in this example
    {NULL, NULL}
} ;

