//------------------------------------------------------------------------------
// LAGraph_A_to_E: Given an undirected graph with no self-loops, builds its
// corresponding incidence matrix
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Vidith Madhu, Texas A&M University

//--------

#include "LG_internal.h"
#include "LAGraphX.h"

// #define dbg

int LAGraph_A_to_E
(
    GrB_Matrix *result, // incidence
    LAGraph_Graph G, // must be undirected, no self-loops
    char *msg
)
{
    // TODO: What are the proper error codes?
    LG_ASSERT_MSG (G->kind == LAGraph_ADJACENCY_UNDIRECTED, -107, "G must be undirected") ;
    LG_ASSERT_MSG (G->nself_edges == 0, -107, "G->nself_edges must be zero") ;

    const GrB_Matrix A = G->A ;
    // Setup
    GrB_Matrix E = NULL ;
    GrB_Index *row_indices = NULL ;
    GrB_Index *col_indices = NULL ;
    // cast everything into a fp64
    double *values = NULL ;

    GrB_Index *E_row_indices = NULL ;
    GrB_Index *E_col_indices = NULL ;
    double *E_values = NULL ;

    GrB_Index nvals ;
    GrB_Index num_nodes ;

    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    GRB_TRY (GrB_Matrix_nrows (&num_nodes, A)) ;
    GRB_TRY (GrB_Matrix_new (&E, GrB_FP64, num_nodes, nvals / 2)) ;

    // Arrays to extract A into
    LG_TRY (LAGraph_Malloc ((void**)(&row_indices), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&col_indices), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&values), nvals, sizeof(double), msg)) ;

    GRB_TRY (GrB_Matrix_extractTuples (row_indices, col_indices, values, &nvals, A)) ;
    #ifdef dbg
        printf("Printing A values\n");
        for (int i = 0; i < nvals; i++){
            printf("%ld %ld %.5f\n", row_indices[i], col_indices[i], values[i]);
        }
    #endif    
    // number of entries in E should be 2 * n_edges
    // n_edges is nvals / 2 (don't count duplicates). So, number of entries in E is just nvals.
    LG_TRY (LAGraph_Malloc ((void**)(&E_row_indices), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&E_col_indices), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&E_values), nvals, sizeof(double), msg)) ;

    // current index in E_* arrays
    GrB_Index pos = 0;

    for (size_t i = 0; i < nvals; i++) {
        // only consider edges if row < col (prevent duplicates)
        GrB_Index row = row_indices[i] ;
        GrB_Index col = col_indices[i] ;
        double value = values[i] ;
        if (row < col) {
            // first put only row values (1st endpoint)
            E_col_indices[pos] = pos ;
            E_row_indices[pos] = row ;
            E_values[pos] = value ;
            // printf("DBG: pos = %lld, [%lld, %lld, %lld]\n", pos, E_col_indices[pos], E_row_indices[pos], E_values[pos]);
            // now put the col values (2nd endpoint)
            E_col_indices[pos + (nvals >> 1)] = pos ;
            E_row_indices[pos + (nvals >> 1)] = col ;
            E_values[pos + (nvals >> 1)] = value ;
            pos++ ;
        }   
    }
    #ifdef dbg
        printf("Printing E values\n");
        for (int i = 0; i < nvals; i++){
            printf("%ld %ld %.5f\n", E_row_indices[i], E_col_indices[i], E_values[i]);
        }
    #endif
    GRB_TRY (GrB_Matrix_build (E, E_row_indices, E_col_indices, E_values, nvals, GrB_SECOND_UINT64)) ;

    LAGraph_Free ((void**)(&row_indices), msg) ;
    LAGraph_Free ((void**)(&col_indices), msg) ;
    LAGraph_Free ((void**)(&values), msg) ;

    LAGraph_Free ((void**)(&E_row_indices), msg) ;
    LAGraph_Free ((void**)(&E_col_indices), msg) ;
    LAGraph_Free ((void**)(&E_values), msg) ;

    (*result) = E ;
    return (GrB_SUCCESS) ;
}