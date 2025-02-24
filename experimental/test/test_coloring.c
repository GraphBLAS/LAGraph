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
    // ------------------------------------------------
    // setup
    // ------------------------------------------------

    /* required initialization (found from other test files) */
    LAGraph_Init(msg);
    LAGraph_Random_Init(msg);

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
    // run algorithm independet set
    // ------------------------------------------------

    GxB_set (GxB_BURBLE, false) ;

    int num_colors = 0;
    double time = LAGraph_WallClockTime();    
    LAGraph_coloring_independent_set_optimized(&C, &num_colors, G, msg);
    time = LAGraph_WallClockTime() - time;

    GxB_set (GxB_BURBLE, false) ;

    printf("\nTook %g seconds\n", time);
    printf("Initial Matrix:\n"); LAGraph_Matrix_Print(G->A, LAGraph_SHORT, stdout, msg);
    printf("Final color vector:\n"); LAGraph_Vector_Print(C, LAGraph_SHORT, stdout, msg);

    
    // ------------------------------------------------
    // check if coloring is valid
    // ------------------------------------------------

    OK (LG_check_coloring(G, C, msg));
    printf("Number of Colors: %d\n", num_colors);


    // ------------------------------------------------
    // run algorithm maximal independet set
    // ------------------------------------------------

    GxB_set (GxB_BURBLE, false) ;

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

    GxB_set (GxB_BURBLE, false) ;

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
    LAGraph_Random_Finalize(msg);
}

TEST_LIST =
{
    {"coloring", test_coloring},
    {NULL, NULL}
};
