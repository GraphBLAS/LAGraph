#include <stdio.h>              // ?
#include "LG_internal.h"        // ?
#include "LAGraph_test.h"       // for TEST_CHECK
#include "LAGraphX.h"           // for LAGraph_coloring_independent_set
#include "LG_Xtest.h"           // for LG_check_coloring


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

#define LEN 512
char filename[LEN + 1];

const char* matrix_files[] = {
    "ldbc-undirected-example-unweighted.mtx",
};

void test_coloring(void)
{
#if LAGRAPH_SUITESPARSE
    // ------------------------------------------------
    // setup
    // ------------------------------------------------

    /* required initialization */
    LAGraph_Init(msg);
    LG_Random_Init(msg);

    /* initializing A (matrix) and C (color vector) */
    GrB_Matrix A = NULL;
    GrB_Vector C = NULL;

    /* open matrix market file */
    snprintf(filename, LEN, LG_DATA_DIR "%s", "ldbc-undirected-example-unweighted.mtx");
    FILE *f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL); // A has been moved into G->A

    // ------------------------------------------------
    // run algorithm independent set
    // ------------------------------------------------

//  OK (LG_SET_BURBLE (false)) ;

    int num_colors = 0;
    double time = LAGraph_WallClockTime();    
    LAGraph_coloring_independent_set(&C, &num_colors, G, msg);
    time = LAGraph_WallClockTime() - time;

    OK (LG_SET_BURBLE (false)) ;

    printf("\nTook %g seconds\n", time);
    printf("Initial Matrix:\n"); LAGraph_Matrix_Print(G->A, LAGraph_SHORT, stdout, msg);
    printf("Final color vector:\n"); LAGraph_Vector_Print(C, LAGraph_SHORT, stdout, msg);

    
    // ------------------------------------------------
    // check if coloring is valid
    // ------------------------------------------------

    OK (LG_check_coloring(G, C, msg));
    printf("Number of Colors: %d\n", num_colors);

    // induce no assigned color error
    #if defined ( COVERAGE )
    GrB_free(&C);
    C = NULL;
    LAGraph_Delete(&G, msg);
    snprintf(filename, LEN, LG_DATA_DIR "%s", "ldbc-undirected-example-unweighted.mtx");
    f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL); // A has been moved into G->A

    GrB_Vector C_dummy;
    GrB_Vector_new (&C_dummy, GrB_UINT64, 1);
    int check_coloring_result = LG_check_coloring(G, C, msg);
    TEST_CHECK (check_coloring_result == -1);
    #endif

    // ------------------------------------------------
    // run algorithm independent set with hack
    // ------------------------------------------------

    // hack the random number generator to induce an error condition
    #if defined ( COVERAGE )
    GrB_free(&C);
    C = NULL;
    LAGraph_Delete(&G, msg);
    snprintf(filename, LEN, LG_DATA_DIR "%s", "ldbc-undirected-example-unweighted.mtx");
    f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL); // A has been moved into G->A

    printf ("Hack the random number generator to induce a stall:\n") ;
    random_hack = true ;
    int result = LAGraph_coloring_independent_set(&C, &num_colors, G, msg);
    random_hack = false ;
    printf ("hack msg: %d %s\n", result, msg) ;
    TEST_CHECK (result == LAGRAPH_CONVERGENCE_FAILURE) ;
    #endif

    // ------------------------------------------------
    // run algorithm maximal independent set
    // ------------------------------------------------

    GrB_free(&C);
    C = NULL;
    LAGraph_Delete(&G, msg);
    
    /* open matrix market file */
    snprintf(filename, LEN, LG_DATA_DIR "%s", "ldbc-undirected-example-unweighted.mtx");
    f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL); // A has been moved into G->A

    printf("Initial Matrix:\n"); LAGraph_Matrix_Print(G->A, LAGraph_SHORT, stdout, msg);

    num_colors = 0;
    time = LAGraph_WallClockTime();    
    LAGraph_coloring_MIS(&C, &num_colors, G, msg);
    time = LAGraph_WallClockTime() - time;

    OK (LG_SET_BURBLE (false)) ;

    printf("\nTook %g seconds\n", time);
    
    printf("Final color vector:\n"); LAGraph_Vector_Print(C, LAGraph_SHORT, stdout, msg);

    
    // ------------------------------------------------
    // check if coloring is valid
    // ------------------------------------------------

    OK (LG_check_coloring(G, C, msg));
    printf("Number of Colors: %d\n", num_colors);


    /* clean up (don't understand this) */
    OK(LAGraph_Delete(&G, msg));
    LAGraph_Finalize(msg);
    LG_Random_Finalize(msg);
#endif
}

TEST_LIST =
{
    {"coloring", test_coloring},
    {NULL, NULL}
};
