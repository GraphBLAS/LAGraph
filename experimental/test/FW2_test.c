//------------------------------------------------------------------------------
// FW/test.c: test Floyd-Warshall method: all pairs shortest paths
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <LAGraph.h>
#include <LAGraphX.h>

// Used as global variable, otherwise run out of memory when creating second graph
// of size VxV to keep track of floyd warshall distances.
int **graph;//[argv[3]][argv[3]];

void floydWarshall(int V) {
    for (int k = 0; k < V; k++) {
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (graph[i][k] == INT_MAX || graph[k][j] == INT_MAX)
                    continue;
                if (graph[i][j] > graph[i][k] + graph[k][j])
                    graph[i][j] = graph[i][k] + graph[k][j];
            }
        }
    }
}

int **parents;

void floydWarshallParents(int V) {
    parents = (int **) malloc(V * sizeof(int*));
    for (int i = 0; i < V; i++) {
        parents[i] = (int *)malloc(V * sizeof(int));
        for (int j = 0; j < V; j++) {
            // if (i == j)
            //     parents[i][j] = -1;
            // else if (i != j)
                parents[i][j] = i + 1;
        }
    }

    for (int k = 0; k < V; k++) {
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (graph[i][k] == INT_MAX || graph[k][j] == INT_MAX)
                    continue;
                if (graph[i][j] > graph[i][k] + graph[k][j]) {
                    graph[i][j] = graph[i][k] + graph[k][j];
                    parents[i][j] = parents[k][j];
                }
            }
        }
    }
}

#define LAGraph_FREE_ALL                        \
    GrB_free (&A) ;                             \
    GrB_free (&Output) ;                        \
    GrB_free (&regResult) ;                     \
    if (graph != NULL)                          \
    {                                           \
        for (int i = 0; i < V; i++)             \
        {                                       \
            free(graph[i]);                     \
        }                                       \
        free(graph);                            \
    }

//****************************************************************************
int main(int argc, char *argv[])
{
    GrB_Matrix A = NULL, Output = NULL ;
    GrB_Type A_type = NULL, Output_type = NULL;
    GrB_Matrix regResult = NULL;
    GrB_Info info ;
    graph = NULL;

    if (argc < 3)
        return -1;

    LAGraph_Init(NULL);

    int V = atoi(argv[2]);
    graph = (int **) malloc(V * sizeof(int*));
    for (int i = 0; i < V; i++)
        graph[i] = (int *)malloc(V * sizeof(int));

    LAGRAPH_OK (GrB_Matrix_new(&A, GrB_FP64, V, V));

    FILE *file;
    file = fopen(argv[1], "r");

    if (file == NULL) {
        printf("file not open\n");
        return 0;
    }

    LAGRAPH_OK (LAGraph_MMRead(&A, &A_type, file, NULL));
    GrB_Index nrows, ncols;

    LAGRAPH_OK (GrB_Matrix_nrows(&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols(&ncols, A));

    int val;
    for (GrB_Index i = 0; i < V; i++) {
        for (GrB_Index j = 0; j < V; j++) {
            LAGRAPH_OK (GrB_Matrix_extractElement(&val, A, i, j));
            if (info == GrB_SUCCESS)
                graph[i][j] = val;
            else
                graph[i][j] = INT_MAX;
        }
    }

    double tic[2], t1, t2 ;

    LAGraph_Tic (tic, NULL);
    floydWarshall(V);
    LAGraph_Toc (&t1, tic, NULL);
    printf("Non-GraphBLAS Floyd Warshall time in seconds: %14.6f\n", t1);

    LAGraph_Tic (tic, NULL);
    LAGraph_FW(A, &Output, &Output_type);
    LAGraph_Toc (&t2, tic, NULL);
    printf("GraphBLAS Floyd Warshall time in seconds:     %14.6f\n", t2);

    LAGRAPH_OK (GrB_Matrix_new(&regResult, Output_type, V, V));
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (graph[i][j] != INT_MAX)
                GrB_Matrix_setElement(regResult, graph[i][j], i, j);
        }
    }

    bool isSame;
    LAGRAPH_OK(LAGraph_IsEqual_type(&isSame, regResult,
                                    Output, Output_type, NULL));
    if (isSame)
        printf("Test passed for file: %s\n\n", argv[1]);
    else
        printf("Test failed for file: %s\n\n", argv[1]);

    LAGraph_FREE_ALL;
    LAGraph_Finalize(NULL);
    return 0;
}
