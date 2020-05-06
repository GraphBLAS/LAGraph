//------------------------------------------------------------------------------
// LAGraph_dense_relabel: dense relabeling of ids to matrix indices
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

// LAGraph_dense_relabel: relabel sparse IDs to dense row/column indices
// Contributed by Marton Elekes and Gabor Szarnyas,
// Budapest University of Technology and Economics
// (with accented characters: M\'{a}rton Elekes and G\'{a}bor Sz\'{a}rnyas.

// Converts array of sparse IDs (ids) to row/column indices between 0...(nids-1).
//
// Gives back two binary matrices for conversion between ID- and index-based vertices.
// id2index vector can be used to look up for indices of chosen IDs.
// id_dimension gives back the height of Id2index matrix and id2index vector. (Same as width of Index2id_handle matrix.)
//   id_dimension is the size that can store the largest ID in the array.
//   Currently it is the largest valid dimension in SuiteSparse:GraphBLAS (GB_INDEX_MAX = 2^60)
//
// Find usage example in /Test/DenseRelabel/denserelabeltest.c

#include "LAGraph_internal.h"
#include <string.h>

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (Id2index_handle) ;    \
    GrB_free (Index2id_handle) ;    \
    GrB_free (id2index_handle) ;    \
    LAGRAPH_FREE (indices) ;        \
    LAGRAPH_FREE (true_values) ;    \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_dense_relabel               // relabel sparse IDs to dense row/column indices
        (
                GrB_Matrix *Id2index_handle, // output matrix: A(id, index)=1 (unfilled if NULL)
                GrB_Matrix *Index2id_handle, // output matrix: B(index, id)=1 (unfilled if NULL)
                GrB_Vector *id2index_handle, // output vector: v(id)=index (unfilled if NULL)
                const GrB_Index *ids,        // array of unique identifiers (under GB_INDEX_MAX=2^60)
                GrB_Index nids,              // number of identifiers
                GrB_Index *id_dimension      // number of rows in Id2index matrix, id2index vector (unfilled if NULL)
        ) {

    GrB_Index *indices = NULL;
    bool *true_values = NULL;

    // from LAGraph_1_to_n.c
    int nthreads = LAGraph_get_nthreads();
    nthreads = LAGRAPH_MIN ((int64_t) (nids / 4096), nthreads);
    nthreads = LAGRAPH_MAX (nthreads, 1);

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!Id2index_handle && !Index2id_handle && !id2index_handle) {
        LAGRAPH_ERROR ("All output mapping arguments are NULL", GrB_NULL_POINTER);
    }
    if (!ids) {
        LAGRAPH_ERROR ("ids is NULL", GrB_NULL_POINTER);
    }

    // the largest valid dimension in SuiteSparse:GraphBLAS (GB_INDEX_MAX)
    GrB_Index id_max_dimension = ((GrB_Index) (1ULL << 60));
    if (id_dimension)
        *id_dimension = id_max_dimension;

    // set indices 0..(nids-1)
    indices = LAGraph_malloc(nids, sizeof(*indices));
    if (!indices) {
        LAGRAPH_ERROR ("Out of Memory", GrB_OUT_OF_MEMORY);
    }
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (size_t i = 0; i < nids; ++i) {
        indices[i] = i;
    }

    // build vector id2index(original_id) = index
    if (id2index_handle) {
        LAGr_Vector_new(id2index_handle, GrB_UINT64, id_max_dimension);
        LAGr_Vector_build(*id2index_handle, ids, indices, nids, GrB_SECOND_UINT64);
    }

    if (Id2index_handle || Index2id_handle) {
        // initialize true values of the matrix
        true_values = LAGraph_malloc(nids, sizeof(*true_values));
        if (!true_values) {
            LAGRAPH_ERROR ("Out of Memory", GrB_OUT_OF_MEMORY);
        }
        memset(true_values, true, nids * sizeof(*true_values));

        // build matrix Index2id(index, original_id) = 1
        if (Index2id_handle) {
            LAGr_Matrix_new(Index2id_handle, GrB_BOOL, nids, id_max_dimension);
            LAGr_Matrix_build(*Index2id_handle, indices, ids, true_values, nids, GrB_SECOND_UINT64);
        }

        // build matrix Id2index(original_id, index) = 1
        if (Id2index_handle) {
            LAGr_Matrix_new(Id2index_handle, GrB_BOOL, id_max_dimension, nids);
            LAGr_Matrix_build(*Id2index_handle, ids, indices, true_values, nids, GrB_SECOND_UINT64);
        }
    }

    return GrB_SUCCESS;
}
