//------------------------------------------------------------------------------
// LAGraph_dense_relabel: dense relabeling of ids to matrix indices
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// LAGraph_dense_relabel: relabel sparse IDs to dense row/column indices
// Contributed by Marton Elekes and Gabor Szarnyas,
// Budapest University of Technology and Economics
// (with accented characters: M\'{a}rton Elekes and G\'{a}bor Sz\'{a}rnyas.

// Converts array of sparse IDs (ids) to row/column indices between 0...(nids-1).
// The order of IDs is kept, therefore ids can be used for index -> ID
// conversion: ids[index]=id.
//
// Gives back two binary matrices for conversion between ID- and index-based
// vertices.
// id2index vector can be used to look up for indices of chosen IDs.
// id_dimension gives back the height of Id2index matrix and id2index vector.
// (Same as width of Index2id_handle matrix.)
//   id_dimension is the size that can store the largest ID in the array.
//   Currently it is the largest valid dimension in GraphBLAS
//  (GrB_INDEX_MAX+1)
//
// Find usage example in /Test/DenseRelabel/dense_relabel_test.c

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_Free ((void **)&indices) ;            \
    LAGraph_Free ((void **)&true_values) ;         \
}

#include <string.h>       // for memset

#include <LAGraph.h>
#include <LAGraphX.h>
#include <LG_internal.h>  // from src/utility

// These should be freed by the calling code
//    GrB_free (Id2index_handle) ;  \
//    GrB_free (Index2id_handle) ;    \
//    GrB_free (id2index_handle) ;    \

//------------------------------------------------------------------------------

GrB_Info LAGraph_dense_relabel   // relabel sparse IDs to dense row/column indices
(
    GrB_Matrix *Id2index_handle, // output matrix: A(id, index)=1 (unfilled if NULL)
    GrB_Matrix *Index2id_handle, // output matrix: B(index, id)=1 (unfilled if NULL)
    GrB_Vector *id2index_handle, // output vector: v(id)=index (unfilled if NULL)
    const GrB_Index *ids,        // array of unique identifiers (<= GrB_INDEX_MAX)
    GrB_Index nids,              // number of identifiers
    GrB_Index *id_dimension      // number of rows in Id2index matrix, id2index vector (unfilled if NULL)
)
{
    GrB_Info info; // needed for LAGRAPH_OK
    GrB_Index *indices = NULL;
    bool *true_values = NULL;

    // from LAGraph_1_to_n.c
    int nthreads;
    LAGRAPH_OK (LAGraph_GetNumThreads(&nthreads, NULL)) ;
    nthreads = LAGraph_MIN ((int64_t) (nids / 4096), nthreads);
    nthreads = LAGraph_MAX (nthreads, 1);

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!Id2index_handle && !Index2id_handle && !id2index_handle) {
        LAGRAPH_ERROR ("All output mapping arguments are NULL", GrB_NULL_POINTER);
    }
    if (!ids) {
        LAGRAPH_ERROR ("ids is NULL", GrB_NULL_POINTER);
    }

    // the largest valid dimension in GraphBLAS
    GrB_Index id_max_dimension = GrB_INDEX_MAX+1;
    if (id_dimension)
        *id_dimension = id_max_dimension;

    // set indices 0..(nids-1)
    indices = LAGraph_Malloc(nids, sizeof(*indices));
    if (!indices) {
        LAGRAPH_ERROR ("Out of Memory", GrB_OUT_OF_MEMORY);
    }

#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (size_t i = 0; i < nids; ++i) {
        indices[i] = i;
    }

    // build vector id2index(original_id) = index
    if (id2index_handle) {
        LAGRAPH_OK (GrB_Vector_new(id2index_handle, GrB_UINT64,
                                   id_max_dimension));
        LAGRAPH_OK (GrB_Vector_build(*id2index_handle,
                                     ids, indices, nids,
                                     GrB_SECOND_UINT64));
    }

    if (Id2index_handle || Index2id_handle) {
        // initialize true values of the matrix
        true_values = LAGraph_Malloc(nids, sizeof(*true_values));
        if (!true_values) {
            LAGRAPH_ERROR ("Out of Memory", GrB_OUT_OF_MEMORY);
        }
        memset(true_values, true, nids * sizeof(*true_values));

        // build matrix Index2id(index, original_id) = 1
        if (Index2id_handle) {
            LAGRAPH_OK (GrB_Matrix_new(Index2id_handle, GrB_BOOL,
                                       nids, id_max_dimension));
            LAGRAPH_OK (GrB_Matrix_build(*Index2id_handle,
                                         indices, ids, true_values, nids,
                                         GrB_SECOND_UINT64));
        }

        // build matrix Id2index(original_id, index) = 1
        if (Id2index_handle) {
            LAGRAPH_OK (GrB_Matrix_new(Id2index_handle, GrB_BOOL,
                                       id_max_dimension, nids));
            LAGRAPH_OK (GrB_Matrix_build(*Id2index_handle,
                                         ids, indices, true_values, nids,
                                         GrB_SECOND_UINT64));
        }
    }

    LAGraph_FREE_ALL;
    return GrB_SUCCESS;
}
