
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

// Contributed by Scott Kolodziej, Texas A&M University

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

GrB_Info LAGraph_sssp11a         // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in
                               // semiring (INT32)
    GrB_Index source,          // source vertex from which to compute
                               // shortest paths
    int32_t delta,             // delta value for delta stepping
    bool AIsAllPositive        // A boolean indicating whether the entries of
                               // matrix A are all positive
);
GrB_Info LAGraph_sssp11b         // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in
                               // semiring (INT32)
    GrB_Index source,          // source vertex from which to compute
                               // shortest paths
    int32_t delta,             // delta value for delta stepping
    bool AIsAllPositive        // A boolean indicating whether the entries of
                               // matrix A are all positive
);


GrB_Info LAGraph_sssp2         // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in semiring (INT32)
    GrB_Index source,          // source vertex from which to compute shortest paths
    int32_t delta               // delta value for delta stepping
);

GrB_Info LAGraph_sssp12a        // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in
                               // semiring (INT32)
    GrB_Index source,          // source vertex from which to compute
                               // shortest paths
    int32_t delta,             // delta value for delta stepping

    // TODO: make this an enum:
    //      case 0: A can have negative, zero, or positive entries
    //      case 1: A can have zero or positive entries
    //      case 2: A only has positive entries (see FIXME below)
    bool AIsAllPositive        // A boolean indicating whether the entries of
                               // matrix A are all positive
);


GrB_Info LAGraph_sssp12b        // single source shortest paths
(
    GrB_Vector *path_length,   // path_length(i) is the length of the shortest
                               // path from the source vertex to vertex i
    GrB_Matrix A,              // input graph, treated as if boolean in
                               // semiring (INT32)
    GrB_Index source,          // source vertex from which to compute
                               // shortest paths
    int32_t delta,             // delta value for delta stepping

    // TODO: make this an enum:
    //      case 0: A can have negative, zero, or positive entries
    //      case 1: A can have zero or positive entries
    //      case 2: A only has positive entries (see FIXME below)
    bool AIsAllPositive        // A boolean indicating whether the entries of
                               // matrix A are all positive
);

