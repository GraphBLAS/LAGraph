//------------------------------------------------------------------------------
// LAGraph/Test/DenseRelabel/dense_relabel_test.c: test program for LAGraph_dense_relabel
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Contributed by Marton Elekes and Gabor Szarnyas, BME

// Usage:
//
// dense_relabel_test

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL            \
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
        LAGraph_FREE_ALL;                                   \
        exit(EXIT_FAILURE);                                 \
    }                                                       \
}


int main(void)
{
    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------
    GrB_Info info; // needed for LAGRAPH_OK
    GrB_Matrix Id2index = NULL;
    GrB_Matrix Index2id = NULL;
    GrB_Vector id2index = NULL;
    GrB_Vector id_vec = NULL;
    GrB_Vector index_vec = NULL;
    GrB_Vector ref_id_vec = NULL;
    GrB_Vector ref_index_vec = NULL;

    LAGRAPH_OK (LAGraph_Init (NULL)) ;

    // prepare array of IDs
    const GrB_Index big_id = 1ULL << 48;
    const GrB_Index index_of_big_id = 2;
    const GrB_Index identifiers[] = {42, 0, big_id, 1};
    const GrB_Index nids = sizeof(identifiers) / sizeof(identifiers[0]);

    // build mappings
    GrB_Index id_dimension;
    LAGRAPH_OK(LAGraph_dense_relabel(&Id2index, &Index2id, &id2index,
                                     identifiers, nids, &id_dimension));

#if !defined(NDEBUG) && defined(LG_SUITESPARSE)
    LAGRAPH_OK(GxB_fprint(Id2index, GxB_COMPLETE, stdout));
    LAGRAPH_OK(GxB_fprint(Index2id, GxB_COMPLETE, stdout));
    LAGRAPH_OK(GxB_fprint(id2index, GxB_COMPLETE, stdout));
#endif

    //--------------------------------------------------------------------------
    // use id2index vector (original_id -> index)
    //--------------------------------------------------------------------------
    GrB_Index index = 0;
    LAGRAPH_OK (GrB_Vector_extractElement(&index, id2index, big_id));
    ASSERT_TRUE(index_of_big_id == index);

    //--------------------------------------------------------------------------
    // use Id2index (original_id -> index)
    //--------------------------------------------------------------------------
    LAGRAPH_OK (GrB_Vector_new(&id_vec, GrB_BOOL, id_dimension));
    LAGRAPH_OK (GrB_Vector_setElement(id_vec, true, big_id));
#if !defined(NDEBUG) && defined(LG_SUITESPARSE)
    LAGRAPH_OK(GxB_fprint(id_vec, GxB_COMPLETE, stdout));
#endif

    LAGRAPH_OK (GrB_Vector_new(&index_vec, GrB_BOOL, nids));
    LAGRAPH_OK (GrB_vxm(index_vec, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL,
                        id_vec, Id2index, GrB_NULL));
#if !defined(NDEBUG) && defined(LG_SUITESPARSE)
    LAGRAPH_OK(GxB_fprint(index_vec, GxB_COMPLETE, stdout));
#endif

    // test
    LAGRAPH_OK (GrB_Vector_new(&ref_index_vec, GrB_BOOL, nids));
    LAGRAPH_OK (GrB_Vector_setElement(ref_index_vec, true, index_of_big_id));
    {
        bool isequal = false;
        LAGRAPH_OK(LAGraph_Vector_IsEqual_type(&isequal, index_vec, ref_index_vec,
                                               GrB_BOOL, NULL));
        ASSERT_TRUE(isequal);
    }

    //--------------------------------------------------------------------------
    // use Index2id (index -> original_id)
    //--------------------------------------------------------------------------
    LAGRAPH_OK (GrB_Vector_clear(id_vec));
    LAGRAPH_OK (GrB_vxm(id_vec, GrB_NULL, GrB_NULL, GxB_LOR_LAND_BOOL, index_vec, Index2id, GrB_NULL));
#if !defined(NDEBUG) && defined(LG_SUITESPARSE)
    LAGRAPH_OK(GxB_fprint(id_vec, GxB_COMPLETE, stdout));
#endif

    // test
    LAGRAPH_OK (GrB_Vector_new(&ref_id_vec, GrB_BOOL, id_dimension));
    LAGRAPH_OK (GrB_Vector_setElement(ref_id_vec, true, big_id));
    {
        bool isequal = false;
        LAGRAPH_OK(LAGraph_Vector_IsEqual_type(&isequal, id_vec, ref_id_vec,
                                               GrB_BOOL, NULL));
        ASSERT_TRUE(isequal);
    }

    LAGraph_FREE_ALL;
    LAGraph_Finalize(NULL);
    return 0;
}
