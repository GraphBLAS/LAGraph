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
    GrB_free(&parent_cpy) ;                                 \
    GrB_free(&parent_sorted) ;                              \
    GrB_free(&sorted_permutation) ;                         \
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
    GrB_Vector parent_cpy = NULL ;          // don't modify input parent
    GrB_Vector parent_sorted = NULL ;
    GrB_Vector sorted_permutation = NULL ;
    GrB_Matrix S = NULL ;
    
    GrB_Index nvals ;

    GrB_Index n ;
    GRB_TRY (GrB_Vector_nvals (&n, parent)) ;

    GRB_TRY (GrB_Vector_dup (&parent_cpy, parent)) ;

    if (!preserve_mapping) {
        // result dim: n' by n (n' is the new number of vertices)
        GrB_Index n_new = 0 ;
        // n_new (or n') is number of unique elements in parent
        // equals grb_reduce with an "equal" operator

        // important step: relabel the vertices from [0...(n - 1)] to [0...(n' - 1)]
        // need to sort to ensure ordering of new labels matches ordering of old labels
        GRB_TRY (GrB_Vector_new (&parent_sorted, GrB_UINT64, n)) ;
        GRB_TRY (GrB_Vector_new (&sorted_permutation, GrB_UINT64, n)) ;
        GRB_TRY (GxB_Vector_sort (parent_sorted, sorted_permutation, GrB_LT_UINT64, parent_cpy, NULL)) ;
        
        GrB_Index prv = UINT64_MAX ;

        // perform relabelling
        for (GrB_Index i = 0; i < n; i++) {
            uint64_t val, perm_idx ;
            GRB_TRY (GrB_Vector_extractElement (&val, parent_sorted, i)) ;
            GRB_TRY (GrB_Vector_extractElement (&perm_idx, sorted_permutation, i)) ;
            if (i > 0) {
                n_new += (val != prv) ;
            }
            GRB_TRY (GrB_Vector_setElement (parent_cpy, n_new, perm_idx)) ;
            prv = val ;
        }
        // since vertices are 0-indexed
        n_new += (n > 0) ;

        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n_new, n)) ;

        GRB_TRY (GrB_free (&parent_sorted)) ;
        GRB_TRY (GrB_free (&sorted_permutation)) ;
    } else {
        // result dim: n by n
        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n, n)) ;
    }

    GrB_Index J = NULL, I = NULL, J_size, I_size, nvals ;
    bool iso ;
    GxB_Vector_unpack_CSC (parent_cpy, &J, &I, &J_size, &I_size, NULL,
        &nvals, NULL) ;
    GrB_Scalar one ;
    GrB_Scalar_new (&one, GrB_UINT64) ;
    GrB_Scalar_setElement (one, 1) ;
    GxB_Matrix_build_Scalar (S, I, J, one, nvals, NULL) ;

    LAGraph_Free ((void **) &I, NULL) ;
    LAGraph_Free ((void **) &J, NULL) ;
    GrB_free (&one) ;

#if 0
GrB_Info GxB_Vector_unpack_CSC  // unpack a CSC vector
(
    GrB_Vector v,       // vector to unpack (type and length unchanged)
    GrB_Index **vi,     // indices
    void **vx,          // values
    GrB_Index *vi_size, // size of vi in bytes
    GrB_Index *vx_size, // size of vx in bytes
    bool *iso,          // if true, v is iso
    GrB_Index *nvals,   // # of entries in vector
    bool *jumbled,      // if true, indices may be unsorted
    const GrB_Descriptor desc
) ;

GrB_Info GxB_Matrix_build_Scalar    // build a matrix from (I,J,scalar) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    GrB_Scalar scalar,              // value for all tuples
    GrB_Index nvals                 // number of tuples
) ;

GrB_Info GrB_Matrix_build_UINT64    // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const uint64_t *X,              // array of values of tuples
    GrB_Index nvals,                // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;
#endif


#if 0
    for (GrB_Index idx = 0; idx < n; idx++) {
        uint64_t i ;
        GRB_TRY (GrB_Vector_extractElement (&i, parent_cpy, idx)) ;
        GrB_Index j = idx ;
        // printf("set (%lld, %lld) to 1\n", i, j);    
        GRB_TRY (GrB_Matrix_setElement (S, 1, i, j)) ;
    }
#endif

    (*result) = S ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
