#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <LAGraph.h>
#include <LG_internal.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;

#define LEN 512
char filename [LEN+1] ;

void test_argminmax (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL;
    GrB_Matrix x = NULL, p = NULL;

    int dim  = 0;
    bool is_min =1;

    //--------------------------------------------------------------------------
    // test with the A matrix
    //--------------------------------------------------------------------------

    // create the graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "structure.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    printf ("\nInput of Matrix:\n") ;
    GxB_print(A, 2);
    // test the algorithm
    OK (LAGraph_argminmax (&x,&p, A,dim,is_min, msg));
    printf("\n") ;
    GxB_print(x,3);
    GxB_print(p,3);
    // print the result
    

    // OK (LAGraph_Matrix_Print (A, LAGraph_COMPLETE, stdout, msg)) ;

    // check the result (ensure Y is equal to G->A)
    // bool ok ;
    // OK (LAGraph_Matrix_IsEqual (&ok, Y, G->A, msg)) ;
    // TEST_CHECK (ok) ;

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    OK (GrB_free (&x)) ;
    OK (GrB_free (&p)) ;

    OK (LAGraph_Delete (&G, msg)) ;

    LAGraph_Finalize (msg) ;
}

//----------------------------------------------------------------------------
// the make program is created by acutest, and it runs a list of tests:
//----------------------------------------------------------------------------

TEST_LIST =
{
    {"Argminmax", test_argminmax},    // just one test in this example
    {NULL, NULL}
} ;
