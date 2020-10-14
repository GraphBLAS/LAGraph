//------------------------------------------------------------------------------
// FW/main.c: test Floyd-Warshall method: all pairs shortest paths
//------------------------------------------------------------------------------

/*

    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.
*/

//------------------------------------------------------------------------------

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
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
