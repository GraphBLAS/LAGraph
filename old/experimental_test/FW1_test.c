//------------------------------------------------------------------------------
// FW/main.c: test Floyd-Warshall method: all pairs shortest paths
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL        \
    GrB_free (&A);              \
    GrB_free (&Output) ;

//****************************************************************************
int main(int argc, char *argv[]) {
    GrB_Matrix A = NULL, Output = NULL;
    GrB_Type A_type = NULL, Output_type = NULL;
    GrB_Info info;

    if (argc < 4)
        return -1;

    LAGraph_Init(NULL);
    int V = atoi(argv[3]);
    GrB_Matrix_new(&A, GrB_FP32, V, V);
    FILE *file;
    file = fopen(argv[1], "r");

    if (file == NULL) {
        printf("file not open\n");
        return 0;
    }

    LAGraph_MMRead(&A, &A_type, file, NULL);

    double tic [2], t;
    LAGraph_Tic(tic, NULL);
    LAGraph_FW(A, &Output, &Output_type);
    LAGraph_Toc(&t, tic, NULL);
    printf("GraphBLAS Floyd Warshall time in seconds: %14.6f\n", t);

    FILE *outputFile;
    outputFile = fopen(argv[2], "w");
    LAGraph_MMWrite_type(Output, Output_type, outputFile, NULL, NULL);

    LAGraph_FREE_ALL;
    LAGraph_Finalize(NULL);

    return 0;
}
