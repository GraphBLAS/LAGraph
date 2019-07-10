#include <stdio.h>
#include <limits.h>
#include <float.h>
// #include "LAGraph.h"
#include "LAGraph_FW.c"
 
#define OK(method) ({                                                                           \
    info = method;                                                                              \
    if (! (info == GrB_SUCCESS || info == GrB_NO_VALUE))                                        \
    {                                                                                           \
        printf ("Error! File: %s line %d [%d] %s\n", __FILE__, __LINE__, info, GrB_error());    \
        FREE_ALL;                                                                               \
        return (info);                                                                          \
    }                                                                                           \
})
 
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
 
#undef FREE_ALL
#define FREE_ALL                                                                            \
    GrB_free (&A) ;                                                                         \
    GrB_free (&Output) ;
 
int main(int argc, char *argv[]) {
    LAGraph_init();
 
    GrB_Matrix A = NULL, Output = NULL ;
    GrB_Info info ;
    int V = atoi(argv[2]);
    graph = (int **) malloc(V * sizeof(int*));
    for (int i = 0; i < V; i++)
        graph[i] = (int *)malloc(V * sizeof(int));
 
    OK (GrB_Matrix_new(&A, GrB_FP64, V, V));
 
    FILE *file;
    file = fopen(argv[1], "r");
 
    if (file == NULL) {
        printf("file not open\n");
        return 0;
    }
 
    OK (LAGraph_mmread(&A, file));
    GrB_Index nrows, ncols;
 
    OK (GrB_Matrix_nrows(&nrows, A));
    OK (GrB_Matrix_ncols(&ncols, A));
 
    int val;
    for (GrB_Index i = 0; i < V; i++) {
        for (GrB_Index j = 0; j < V; j++) {
            OK (GrB_Matrix_extractElement(&val, A, i, j));
            if (info == GrB_SUCCESS)
                graph[i][j] = val;
            else
                graph[i][j] = INT_MAX;
        }
    }
 
    double tic[2], t1, t2 ;
 
    LAGraph_tic (tic);
    floydWarshall(V);
    t1 = LAGraph_toc (tic);
    printf("Non-GraphBLAS Floyd Warshall time in seconds: %14.6f\n", t1);
 
    LAGraph_tic (tic);
    LAGraph_FW(A, &Output);
    t2 = LAGraph_toc (tic);
 
    printf("GraphBLAS Floyd Warshall time in seconds:     %14.6f\n", t2);
 
    GrB_Matrix regResult = NULL;
    OK (GrB_Matrix_new(&regResult, GrB_FP64, V, V));
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (graph[i][j] != INT_MAX)
                GrB_Matrix_setElement(regResult, graph[i][j], i, j);
        }
    }
 
    bool isSame;
    OK(LAGraph_isequal(&isSame, regResult, Output, GrB_EQ_INT32));
    if (isSame)
        printf("Test passed for file: %s\n\n", argv[1]);
    else
        printf("Test failed for file: %s\n\n", argv[1]);
     
    LAGraph_finalize();
    return 0;
}
