//------------------------------------------------------------------------------
// LAGraph/Test/DenseRelabel/dense_relabel_test.c: test program for LAGraph_dense_relabel
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

// Contributed by Marton Elekes and Gabor Szarnyas, BME

// Usage:
//
// dense_relabel_test

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&Id2index) ;          \
    GrB_free (&Index2id) ;          \
    GrB_free (&id2index) ;          \
    GrB_free (&id_vec) ;            \
    GrB_free (&index_vec) ;         \
    GrB_free (&ref_id_vec) ;        \
    GrB_free (&ref_index_vec) ;     \
}

#define ASSERT_TRUE(expr)                                   \
{                                                           \
    if(!(expr)) {                                           \
        fprintf(stderr, "Test failed: %s\nFile: %s:%d\n",   \
                #expr, __FILE__, __LINE__);                 \
        LAGRAPH_FREE_ALL;                                   \
        exit(EXIT_FAILURE);                                 \
    }                                                       \
}


int main(void) {

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Matrix Id2index = NULL;
    GrB_Matrix Index2id = NULL;
    GrB_Vector id2index = NULL;
    GrB_Vector id_vec = NULL;
    GrB_Vector index_vec = NULL;
    GrB_Vector ref_id_vec = NULL;
    GrB_Vector ref_index_vec = NULL;

    LAGRAPH_TRY_CATCH (LAGraph_init());

    // prepare array of IDs
    const GrB_Index big_id = 1ULL << 48;
    const GrB_Index index_of_big_id = 2;
    const GrB_Index identifiers[] = {42, 0, big_id, 1};
    const GrB_Index nids = sizeof(identifiers) / sizeof(identifiers[0]);

    // build mappings
    GrB_Index id_dimension;
    LAGRAPH_TRY_CATCH(LAGraph_dense_relabel(&Id2index, &Index2id, &id2index, identifiers, nids, &id_dimension));

#ifndef NDEBUG
    LAGRAPH_TRY_CATCH(GxB_fprint(Id2index, GxB_COMPLETE, stdout));
    LAGRAPH_TRY_CATCH(GxB_fprint(Index2id, GxB_COMPLETE, stdout));
    LAGRAPH_TRY_CATCH(GxB_fprint(id2index, GxB_COMPLETE, stdout));
#endif

    //--------------------------------------------------------------------------
    // use id2index vector (original_id -> index)
    //--------------------------------------------------------------------------
    GrB_Index index = 0;
    LAGr_Vector_extractElement(&index, id2index, big_id);
    ASSERT_TRUE(index_of_big_id == index);

    //--------------------------------------------------------------------------
    // use Id2index (original_id -> index)
    //--------------------------------------------------------------------------
    LAGr_Vector_new(&id_vec, GrB_BOOL, id_dimension);
    LAGr_Vector_setElement(id_vec, true, big_id);
#ifndef NDEBUG
    LAGRAPH_TRY_CATCH(GxB_fprint(id_vec, GxB_COMPLETE, stdout));
#endif

    LAGr_Vector_new(&index_vec, GrB_BOOL, nids);
    LAGr_vxm(index_vec, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL, id_vec, Id2index, GrB_NULL);
#ifndef NDEBUG
    LAGRAPH_TRY_CATCH(GxB_fprint(index_vec, GxB_COMPLETE, stdout));
#endif

    // test
    LAGr_Vector_new(&ref_index_vec, GrB_BOOL, nids);
    LAGr_Vector_setElement(ref_index_vec, true, index_of_big_id);
    {
        bool isequal = false;
        LAGRAPH_TRY_CATCH(LAGraph_Vector_isequal(&isequal, index_vec, ref_index_vec, GrB_NULL));
        ASSERT_TRUE(isequal);
    }

    //--------------------------------------------------------------------------
    // use Index2id (index -> original_id)
    //--------------------------------------------------------------------------
    LAGr_Vector_clear(id_vec);
    LAGr_vxm(id_vec, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL, index_vec, Index2id, GrB_NULL);
#ifndef NDEBUG
    LAGRAPH_TRY_CATCH(GxB_fprint(id_vec, GxB_COMPLETE, stdout));
#endif

    // test
    LAGr_Vector_new(&ref_id_vec, GrB_BOOL, id_dimension);
    LAGr_Vector_setElement(ref_id_vec, true, big_id);
    {
        bool isequal = false;
        LAGRAPH_TRY_CATCH(LAGraph_Vector_isequal(&isequal, id_vec, ref_id_vec, GrB_NULL));
        ASSERT_TRUE(isequal);
    }

    LAGRAPH_FREE_ALL;
    LAGraph_finalize();
}
