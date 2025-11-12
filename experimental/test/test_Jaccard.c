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


void test_Jaccard (void)
{
    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, ABool = NULL;
	GrB_Index n; 

    // create the graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "zenios.mtx") ;
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
    //FILE *fo = fopen ("LAGraph_Jaccard_all_pairs.txt", "w") ;
//	OK (GxB_Matrix_fprint(JC, "my matrix", LAGraph_COMPLETE, fo));
//	fclose(fo);
	OK(GrB_free(&JC)); 
	JC = NULL;

	//--------------------------------------------------------------------------
    // Experiment 2: Jaccard weights (neighbors only)
    //--------------------------------------------------------------------------
	printf("Exp 2 : Jaccard weights\n");
	all_pairs = false;	// if false it only computes jaccard weights(jaccard similarity for neighbors)
    OK (LAGraph_Jaccard(&JC, G, all_pairs, msg)) ;
//    fo = fopen ("LAGraph_Jaccard_weights.txt", "w") ;
//    OK (GxB_Matrix_fprint(JC, "my matrix", LAGraph_COMPLETE, fo));
//    fclose(fo);
    TEST_CHECK(0 == GrB_free(&JC));

    // free everything
	OK(LAGraph_Delete(&G, msg));
    OK(GrB_free(&A));
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

