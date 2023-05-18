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

#define LG_FREE_WORK                                \
{                                                   \
    GrB_free(&eq) ;                                 \
    GrB_free(&id) ;                                 \
    LAGraph_Free((void**) indices, msg) ;           \
    LAGraph_Free((void**) vals, msg) ;              \
}                                                   \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    GrB_free(&S) ;                                  \
}                                                   \

int LAGraph_Parent_to_S
(
    GrB_Matrix *result,
    GrB_Vector parent,      // dense vector of size n. parent[i] -> representative of node i
    int preserve_mapping,   // whether to preserve the original namespace of nodes
    char *msg
)
{
    GrB_Vector eq = NULL ;
    GrB_Vector id = NULL ;  // id[i] = i
    GrB_Matrix S = NULL ;

    GrB_Index *indices ;
    GrB_Index *vals ;
    GrB_Index nvals ;

    GrB_Index n ;
    GRB_TRY (GrB_Vector_nvals (&n, parent)) ;

    if (!preserve_mapping) {
        // result dim: n' by n (n' != n)
        GrB_Index n_new ;
        // n_new is number of unique elements in parent
        // equals grb_reduce with an "equal" operator
        GRB_TRY (GrB_Vector_new (&id, GrB_UINT64, n)) ;
        GRB_TRY (GrB_Vector_new (&eq, GrB_BOOL, n)) ;

        GRB_TRY (GrB_eWiseAdd (id, NULL, NULL, GxB_FIRSTI_INT64, id, id, NULL)) ;
        GRB_TRY (GrB_eWiseAdd (eq, NULL, NULL, GxB_ISEQ_INT64, parent, id, NULL)) ;

        GRB_TRY (GrB_reduce (&n_new, NULL, GrB_PLUS_MONOID_UINT64, eq, NULL)) ;

        GRB_TRY (GrB_Matrix_new (&S, GrB_BOOL, n_new, n)) ;

        GRB_TRY (GrB_free (&eq)) ;
        GRB_TRY (GrB_free (&id)) ;
    } else {
        // result dim: n by n
        GRB_TRY (GrB_Matrix_new (&S, GrB_BOOL, n, n)) ;
    }

    // extract tuples from parent, make S[j][i] = 1 for p[i] = j

    GRB_TRY (GrB_Vector_nvals (&nvals, parent)) ;

    LG_TRY (LAGraph_Malloc ((void**)(&indices), nvals, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&vals), nvals, sizeof(GrB_Index), msg)) ;

    GRB_TRY (GrB_Vector_extractTuples (indices, vals, &nvals, parent)) ;

    for (GrB_Index idx = 0; idx < nvals; idx++) {
        GrB_Index i = vals[idx] ;
        GrB_Index j = indices[idx] ;
        
        GRB_TRY (GrB_Matrix_setElement (S, true, i, j)) ;
    }

    (*result) = S ;
    return (GrB_SUCCESS) ;
}