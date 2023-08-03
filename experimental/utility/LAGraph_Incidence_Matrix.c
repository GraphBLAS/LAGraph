//------------------------------------------------------------------------------
// LAGraph_Incidence_Matrix: Given the adjacency matrix of an undirected graph with no 
// self-loops, builds its corresponding incidence matrix
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Vidith Madhu, Texas A&M University

// Given an undirected graph G, construct the incidence matrix E.
/*
The incidence matrix E has size n-by-e where the
undirected graph G has n nodes and e edges.  If the kth edge of G is the edge
(i,j), then the column E(:,k) contains two entries:  E(i,k) and E(j,k), which
have the same value.  If the graph G is weighted, then both E(i,k) and E(j,k)
are equal to the weight of the (i,j) edge.  If G is unweighted, then both are
equal to 1 (and the matrix E is thus iso-valued).

The result has the same type as the input adjacency matrix.
*/

//--------

#include "LG_internal.h"
#include "LAGraphX.h"

#include <omp.h>

// #define dbg

#undef LG_FREE_ALL
#define LG_FREE_ALL                                           \
{                                                             \
   LAGraph_Free ((void**)(&row_indices), msg) ;               \
   LAGraph_Free ((void**)(&col_indices), msg) ;               \
   LAGraph_Free ((void**)(&values), msg) ;                    \
   LAGraph_Free ((void**)(&ramp), msg) ;                      \
   GrB_free (&E_half) ;                                       \
   GrB_free (&A_tril) ;                                       \
}                                                             \

int LAGraph_Incidence_Matrix
(
    GrB_Matrix *result, // incidence
    LAGraph_Graph G, // must be undirected, no self-loops
    char *msg
)
{
    
    GrB_Matrix E = NULL ;
    GrB_Matrix E_half = NULL ;
    GrB_Matrix A_tril = NULL ;

    GrB_Index *row_indices = NULL ;
    GrB_Index *col_indices = NULL ;
    void *values = NULL ;

    GrB_Index *ramp = NULL ;

    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED,
        LAGRAPH_INVALID_GRAPH, 
        "G must be undirected"
    ) ;

    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, "G->nself_edges must be zero") ;

    GrB_Matrix A = G->A ;

    char typename[LAGRAPH_MAX_NAME_LEN] ;
    GrB_Type type ;
    LG_TRY (LAGraph_Matrix_TypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;

    GrB_Index nvals ;
    GrB_Index num_nodes, num_edges ;

    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    num_edges = nvals / 2 ;
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;

    GRB_TRY (GrB_Matrix_new (&A_tril, type, num_nodes, num_nodes)) ;
    GRB_TRY (GrB_Matrix_new (&E, type, num_nodes, num_edges)) ;
    GRB_TRY (GrB_Matrix_new (&E_half, type, num_nodes, num_edges)) ;

    bool is_uint64 = (type == GrB_UINT64) ;
    bool is_float = ((type == GrB_FP32) || (type == GrB_FP64)) ;

    // which intermediate type to cast to
    // uint64_t : 0, double: 1, int64_t: 2
    // if the input matrix type is GrB_INT* or GrB_UINT{8|16|32}, it becomes a int64_t
    // if the input matrix type is GrB_UINT64, it becomes a uint64_t
    // if the input matrix type is GrB_FP{32|64}, it becomes a double
    int which_type = (is_uint64 ? 0 : (is_float ? 1 : 2)) ;

    size_t value_sizes[3] = {sizeof(uint64_t), sizeof(double), sizeof(int64_t)} ;

    // (*result) = E ;
    // return (GrB_SUCCESS) ;

    // get just the lower triangular entries
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;

    // Arrays to extract A into
    LG_TRY (LAGraph_Malloc ((void**)(&row_indices), num_edges, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&col_indices), num_edges, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&values), num_edges, value_sizes[which_type], msg)) ;

    switch (which_type){
        case 0:
            GRB_TRY (GrB_Matrix_extractTuples_UINT64 (row_indices, col_indices, values, &num_edges, A_tril)) ;
            break;
        case 1:
            GRB_TRY (GrB_Matrix_extractTuples_FP64 (row_indices, col_indices, values, &num_edges, A_tril)) ;
            break;
        case 2:
            GRB_TRY (GrB_Matrix_extractTuples_INT64 (row_indices, col_indices, values, &num_edges, A_tril)) ;
            break;
    }
    
    #ifdef dbg
        printf("Printing A_tril values\n");
        for (int i = 0; i < num_edges; i++){
            printf("%ld %ld %.5f\n", row_indices[i], col_indices[i], values[i]);
        }
    #endif
    LG_TRY (LAGraph_Malloc ((void**)(&ramp), num_edges, sizeof(GrB_Index), msg)) ;

    #pragma omp parallel for
    for (GrB_Index i = 0 ; i < num_edges ; i++) {
        ramp[i] = i ;
    }

    // build E_1 with (row_indices, ramp, values)
    // build E_2 with (col_indices, ramp, values)
    // E = E_1 + E_2
    switch (which_type) {
        case 0:
            GRB_TRY (GrB_Matrix_build_UINT64 (E_half, col_indices, ramp, values, num_edges, NULL)) ;
            GRB_TRY (GrB_Matrix_build_UINT64 (E, row_indices, ramp, values, num_edges, NULL)) ;
            break;
        case 1:
            GRB_TRY (GrB_Matrix_build_FP64 (E_half, col_indices, ramp, values, num_edges, NULL)) ;
            GRB_TRY (GrB_Matrix_build_FP64 (E, row_indices, ramp, values, num_edges, NULL)) ;
            break;
        case 2:
            GRB_TRY (GrB_Matrix_build_INT64 (E_half, col_indices, ramp, values, num_edges, NULL)) ;
            GRB_TRY (GrB_Matrix_build_INT64 (E, row_indices, ramp, values, num_edges, NULL)) ;
            break;
    }

    GRB_TRY (GrB_eWiseAdd (E, NULL, NULL, GrB_PLUS_FP64, E, E_half, NULL)) ;

    // LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg) ;

    LG_FREE_ALL ;
    
    (*result) = E ;
    return (GrB_SUCCESS) ;
}
