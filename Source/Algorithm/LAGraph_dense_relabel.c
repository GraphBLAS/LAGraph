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

// LAGraph_dense_relabel: Contributed by Marton Elekes and Gabor Szarnyas,
// Budapest University of Technology and Economics
// (with accented characters: M\'{a}rton Elekes and G\'{a}bor Sz\'{a}rnyas.

// ...

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (MMapping_handle) ;    \
    GrB_free (VMapping_handle) ;    \
    LAGRAPH_FREE (indices) ;        \
    LAGRAPH_FREE (true_values) ;    \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_dense_relabel   // compute dense relabel
        (
                GrB_Matrix *MMapping_handle, // output matrix with the mapping (unfilled if NULL)
                GrB_Vector *VMapping_handle, // output vector with the mapping (unfilled if NULL)
                const GrB_Index *ids,        // array of identifiers
                GrB_Index nids               // number of identifiers
        ) {

    GrB_Info info;
    GrB_Index *indices = NULL;
    bool *true_values = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!MMapping_handle && !VMapping_handle) {
        LAGRAPH_ERROR ("Both mapping arguments are NULL", GrB_NULL_POINTER);
    }

    // the largest valid dimension in SuiteSparse:GraphBLAS
    GrB_Index id_max_dimension = ((GrB_Index) (1ULL << 60));

    // set indices 0...nids-1
    indices = LAGraph_malloc(nids, sizeof(*indices));
    if (!indices) {
        LAGRAPH_ERROR ("Out of Memory", GrB_OUT_OF_MEMORY);
    }
    for (size_t i = 0; i < nids; ++i) {
        indices[i] = i;
    }

    // build vector VMapping[original_id] = index
    if (VMapping_handle) {
        LAGRAPH_OK(GrB_Vector_new(VMapping_handle, GrB_UINT64, id_max_dimension));
        LAGRAPH_OK(GrB_Vector_build_UINT64(*VMapping_handle, ids, indices, nids, GrB_SECOND_UINT64));
    }
    if (MMapping_handle) {
        // initilize true values of the matrix
        true_values = LAGraph_malloc(nids, sizeof(*true_values));
        if (!true_values) {
            LAGRAPH_ERROR ("Out of Memory", GrB_OUT_OF_MEMORY);
        }
        for (size_t i = 0; i < nids; ++i) {
            true_values[i] = true;
        }

        // build matrix MMapping[index, original_id] = 1
        LAGRAPH_OK(GrB_Matrix_new(MMapping_handle, GrB_BOOL, nids, id_max_dimension));
        LAGRAPH_OK(GrB_Matrix_build_BOOL(*MMapping_handle, indices, ids, true_values, nids, GrB_SECOND_UINT64));
    }

    return GrB_SUCCESS;
}
