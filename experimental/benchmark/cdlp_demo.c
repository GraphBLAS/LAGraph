#include "LAGraph_demo.h"

#define LG_FREE_ALL                     \
{                                       \
    LAGraph_Delete (&G, msg) ;          \
    GrB_free (&A) ;                     \
    GrB_free (&CDLP_vector) ;           \
}

int main(int argc, char **argv)
{
    char msg[LAGRAPH_MSG_LEN];

    LAGraph_Graph G = NULL;
    GrB_Matrix A = NULL;
    GrB_Vector CDLP_vector = NULL; // CDLP result vector

    // Start GraphBLAS and LAGraph
    bool burble = false;
    demo_init(burble);

    // Read the input matrix
    char *matrix_name = (argc > 1) ? argv[1] : "stdin";
    LAGRAPH_TRY(readproblem(&G, NULL, false, false, false, NULL, true, argc, argv));

    // Prepare parameters for CDLP
    bool symmetric = true; // Set based on your matrix characteristics
    bool sanitize = false; // Set as needed
    int itermax = 100; // Maximum number of iterations
    double t[2]; // To store timings

    // Run the CDLP algorithm
    LAGRAPH_TRY(LAGraph_cdlp(&CDLP_vector, G->A, symmetric, sanitize, itermax, t, msg));

    // Print results and timings
    GxB_print(CDLP_vector, GxB_SHORT);
    printf("Sanitize time: %f\n", t[0]);
    printf("CDLP time: %f\n", t[1]);

    // Optionally, save the result to a file
    // ... [File writing code similar to the original program]

    // Cleanup and finalize
    LG_FREE_ALL;
    LAGRAPH_TRY(LAGraph_Finalize(msg));

    return (GrB_SUCCESS);
}
