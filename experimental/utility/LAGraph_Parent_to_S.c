//------------------------------------------------------------------------------
// LAGraph_Parent_to_S: Given a dense parent vector for an undirected graph, builds the
// corresponding S matrix needed to coarsen the graph
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Vidith Madhu, Texas A&M University

//--------

#include "LG_internal.h"
#include "LAGraphX.h"

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                                        \
{                                                           \
    GrB_free(&parent_sorted) ;                              \
    GrB_free(&sorted_permutation) ;                         \
    LAGraph_Free((void**)(&new_labels), msg) ;              \
    LAGraph_Free ((void**)(&new_label_assigned), msg);      \
}                                                           \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    GrB_free(&S) ;                                  \
}                                                   \

int LAGraph_Parent_to_S
(
    GrB_Matrix *result,
    GrB_Vector parent,      // dense integer vector of size n. parent[i] -> representative of node i
    int preserve_mapping,   // whether to preserve the original namespace of nodes
    char *msg
)
{
    GrB_Vector parent_sorted = NULL ;
    GrB_Vector sorted_permutation = NULL ;
    GrB_Matrix S = NULL ;

    GrB_Index *indices ;
    GrB_Index *vals ;
    GrB_Index *new_labels ;     // tracks new label mapping
    bool *new_label_assigned ;
    GrB_Index nvals ;

    GrB_Index n ;
    GRB_TRY (GrB_Vector_nvals (&n, parent)) ;

    if (!preserve_mapping) {
        // result dim: n' by n (n' is the new number of vertices)
        GrB_Index n_new = 0 ;
        // n_new (or n') is number of unique elements in parent
        // equals grb_reduce with an "equal" operator

        // important step: relabel the vertices from [0...(n - 1)] to [0...(n' - 1)]
        // need to sort to ensure ordering of new labels matches ordering of old labels
        GRB_TRY (GrB_Vector_new (&parent_sorted, GrB_UINT64, n)) ;
        GRB_TRY (GrB_Vector_new (&sorted_permutation, GrB_UINT64, n)) ;
        GRB_TRY (GxB_Vector_sort (parent_sorted, sorted_permutation, GrB_LT_UINT64, parent, NULL)) ;

        LG_TRY (LAGraph_Malloc ((void**)(&new_labels), n, sizeof(GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**)(&new_label_assigned), n, sizeof(bool), msg)) ;
        
        memset (new_label_assigned, 0, sizeof(bool) * n) ;

        // perform relabelling
        for (GrB_Index i = 0; i < n; i++) {
            uint64_t val, perm_idx ;
            GRB_TRY (GrB_Vector_extractElement (&val, parent_sorted, i)) ;
            GRB_TRY (GrB_Vector_extractElement (&perm_idx, sorted_permutation, i)) ;
            if (!new_label_assigned[val]) {
                new_labels[val] = n_new ;
                new_label_assigned[val] = true ;
                n_new++ ;
            }
            GRB_TRY (GrB_Vector_setElement (parent, new_labels[val], perm_idx)) ;
        }

        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n_new, n)) ;

        GRB_TRY (GrB_free (&parent_sorted)) ;
        GRB_TRY (GrB_free (&sorted_permutation)) ;
        LG_TRY (LAGraph_Free ((void**)(&new_labels), msg)) ;
        LG_TRY (LAGraph_Free ((void**)(&new_label_assigned), msg)) ;
    } else {
        // result dim: n by n
        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n, n)) ;
    }

    for (GrB_Index idx = 0; idx < n; idx++) {
        uint64_t i ;
        GRB_TRY (GrB_Vector_extractElement (&i, parent, idx)) ;
        GrB_Index j = idx ;
        // printf("set (%lld, %lld) to 1\n", i, j);    
        GRB_TRY (GrB_Matrix_setElement (S, 1, i, j)) ;
    }

    (*result) = S ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}