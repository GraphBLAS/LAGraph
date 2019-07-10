#include "LAGraph.h"
//#include "simple_timer.h"
#include "LAGraph_FW.c"

 
#undef FREE_ALL
#define FREE_ALL                                                                            \
    GrB_free (&A);
 
int main(int argc, char *argv[]) {
    LAGraph_init();
 
    GrB_init (GrB_NONBLOCKING);
 
    GrB_Matrix A = NULL, Output = NULL;
    GrB_Info info;
    int V = atoi(argv[3]);
    GrB_Matrix_new(&A, GrB_FP32, V, V);
    FILE *file;
    file = fopen(argv[1], "r");
 
    if (file == NULL) {
        printf("file not open\n");
        return 0;
    }
 
    LAGraph_mmread(&A, file);
 
    // double tic [2], t;
    // simple_tic(tic);
    LAGraph_FW(A, &Output);
    // t = simple_toc (tic);
    // printf("GraphBLAS Floyd Warshall time in seconds: %14.6f\n", t);
 
    FILE *outputFile;
    outputFile = fopen(argv[2], "w");
    LAGraph_mmwrite(Output, outputFile);
 
    GrB_finalize();
    LAGraph_finalize();
 
    return 0;
}
