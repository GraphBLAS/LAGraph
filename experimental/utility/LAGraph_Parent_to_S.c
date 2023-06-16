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
    GrB_free(&apply_mask) ;                                 \
    GrB_free (&one) ;                                       \
    GrB_free(&VALUENEQ_ROWINDEX_UINT64) ;                   \
    GrB_free(&ROWINDEX_SUBK_UINT64) ;                       \
    LAGraph_Free ((void**) &discard_indices, NULL) ;        \
    LAGraph_Free ((void**) &discard_values, NULL) ;         \
    LAGraph_Free ((void**) &I, NULL) ;                      \
    LAGraph_Free ((void**) &J, NULL) ;                      \
}                                                           \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    GrB_free(&S) ;                                  \
}                                                   \


#define F_INDEX_UNARY(f)  ((void (*)(void *, const void *, GrB_Index, GrB_Index, const void *)) f)

void valueneq_index_func (bool *z, const uint64_t *x, GrB_Index i, GrB_Index j, const void *y) {
    (*z) = ((*x) != i) ;
}

void index_subk_func (uint64_t *z, const uint64_t *x, GrB_Index i, GrB_Index j, const uint64_t *y) {
    (*z) = (i - (*y)) ;
}

int LAGraph_Parent_to_S
(
    GrB_Matrix *result,
    GrB_Vector parent,      // dense integer vector of size n. parent[i] -> representative of node i
    int preserve_mapping,   // whether to preserve the original namespace of nodes
    char *msg
)
{
    GrB_Vector parent_cpy = NULL ;          // don't modify input parent (also useful to have for compressing node labels)
    GrB_Vector parent_sorted = NULL ;
    GrB_Vector sorted_permutation = NULL ;
    GrB_Matrix S = NULL ;

    // used to unpack parent_cpy and build S
    GrB_Index *J = NULL, *I = NULL, J_size, I_size ;
    GrB_Scalar one = NULL ;

    // used for node compression
    GrB_Index *discard_indices = NULL ;
    uint64_t *discard_values = NULL ;
    GrB_Vector apply_mask = NULL ;
    // GrB_Vector full = NULL ;

    GrB_IndexUnaryOp VALUENEQ_ROWINDEX_UINT64 = NULL ;
    GrB_IndexUnaryOp ROWINDEX_SUBK_UINT64 = NULL ;

    GrB_Index n ;
    GRB_TRY (GrB_Vector_nvals (&n, parent)) ;

    if (preserve_mapping) {
        // parent_cpy will be the same as parent
        GRB_TRY (GrB_Vector_dup (&parent_cpy, parent)) ;
    } else {
        // we want an empty vector for the compression step (see below code)
        GRB_TRY (GrB_Vector_new (&parent_cpy, GrB_UINT64, n)) ;
    }

    if (!preserve_mapping) {
#if 0
        // old code
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
#endif
        /*
        new code:
            - GrB_select to identify discarded nodes
            - suppose we have m discarded nodes - want to do O(m) grb_apply's
            - example:
                parent vector:
                    0 1 2 1 4 7 6 7
                    . . . ^ . ^ . . (discarded nodes marked with '^')
                claim: All other nodes besides those discarded are in order
                proof: If not discarded, then p[i] = i, which is in order w.r.t. index
            
            - So, we can create the new mapping while preserving order without a sort, and
              without looping through all nodes
            - The looping is only O(m) for m discarded nodes.
                - In the case of coarsening, if we sum this across all coarsening steps its just O(n).
        */
        
        GRB_TRY (GrB_IndexUnaryOp_new (&VALUENEQ_ROWINDEX_UINT64, F_INDEX_UNARY(valueneq_index_func), GrB_BOOL, GrB_UINT64, GrB_UINT64)) ;
        GRB_TRY (GrB_IndexUnaryOp_new (&ROWINDEX_SUBK_UINT64, F_INDEX_UNARY(index_subk_func), GrB_UINT64, GrB_UINT64, GrB_UINT64)) ;
        GRB_TRY (GrB_Vector_new (&apply_mask, GrB_UINT64, n)) ;
        // GRB_TRY (GrB_Vector_new (&full, GrB_UINT64, n)) ;
        // GRB_TRY (GrB_assign (full, NULL, NULL, 1, GrB_ALL, n, NULL)) ;

        // put all discarded nodes into parent_cpy
        GRB_TRY (GrB_select (parent_cpy, NULL, NULL, VALUENEQ_ROWINDEX_UINT64, parent, 0, NULL)) ;
        // now, parent_cpy has entries for all discarded nodes
        // we want to get the indices of these entries
        GrB_Index num_discarded ;
        GRB_TRY (GrB_Vector_nvals (&num_discarded, parent_cpy)) ;
        LG_TRY (LAGraph_Malloc ((void**)(&discard_indices), num_discarded, sizeof(GrB_Index), msg)) ;
        LG_TRY (LAGraph_Malloc ((void**)(&discard_values), num_discarded, sizeof(uint64_t), msg)) ;

        GRB_TRY (GrB_Vector_extractTuples (discard_indices, discard_values, &num_discarded, parent_cpy)) ;

        // now, perform our GrB_apply's to do the relabeling
        for (GrB_Index i = 0 ; i < num_discarded ; i++) {
            GrB_Index discarded_node = discard_indices [i] ;
            if (i == 0) {
                // first apply everything
                GRB_TRY (GrB_apply (parent_cpy, NULL, NULL, GrB_IDENTITY_UINT64, parent, NULL)) ;
            }
            // mask everything up to and including discarded_node
            GRB_TRY (GrB_select (apply_mask, NULL, NULL, GrB_ROWLE, parent, discarded_node, NULL)) ;

            // now do the GrB_apply using the complement of apply_mask
            // how this works: suppose we are at the ith non-discarded entry and want to compute its new value.
            // we know that p[i] = i. If there are k discarded entries before the ith entry, then this will remap to p[i] = i - k.
            GRB_TRY (GrB_apply (parent_cpy, apply_mask, NULL, ROWINDEX_SUBK_UINT64, parent, i + 1, GrB_DESC_SC)) ;
        }

        // update the entries for discarded nodes
        // first, store the new entries (don't intermix extractElement and setElement)
        for (GrB_Index i = 0 ; i < num_discarded ; i++) {
            GrB_Index old_parent = discard_values [i] ;
            uint64_t new_parent ;
            GRB_TRY (GrB_Vector_extractElement (&new_parent, parent_cpy, old_parent)) ;
            discard_values [i] = new_parent ;
        }

        for (GrB_Index i = 0 ; i < num_discarded ; i++) {
            GrB_Index discarded_node = discard_indices [i] ;
            GRB_TRY (GrB_Vector_setElement (parent_cpy, discard_values [i], discarded_node)) ;
        }

        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n - num_discarded, n)) ;
    } else { 
        // result dim: n by n
        GRB_TRY (GrB_Matrix_new (&S, GrB_UINT64, n, n)) ;
    }

    GrB_Index nvals ;

    GRB_TRY (GxB_Vector_unpack_CSC (parent_cpy, &J, (void**) &I, &J_size, &I_size, NULL, &nvals, NULL, NULL)) ;
    GRB_TRY (GrB_Scalar_new (&one, GrB_UINT64)) ;
    GRB_TRY (GrB_Scalar_setElement (one, 1)) ;

    GRB_TRY (GxB_Matrix_build_Scalar (S, I, J, one, nvals)) ;    

#if 0
    // old code (not parallel)
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
 