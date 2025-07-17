//------------------------------------------------------------------------------
// LAGraph_Matrix_Hash: generate a single hash value for an entire matrix
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>

#undef LG_FREE_ALL
#define LG_FREE_ALL               \
{                                 \
    GrB_free(&C);                 \
    GrB_free(&lg_hash_edge);      \
}

#define GOLDEN_GAMMA 0x9E3779B97F4A7C15LL

// The init function computes a cheesy hash based on splitmix64t
void LG_HM_hash_edge (uint64_t *z, const uint64_t *x,
    GrB_Index i, GrB_Index j, const uint64_t *seed)
{
    uint64_t result = (i + j * i  + GOLDEN_GAMMA) ;
    result = (result ^ (*x)) * 0xBF58476D1CE4E5B9LL ;
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9LL ;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EBLL ;
    result = (result ^ (result >> 31)) ;
    (*z) = result ;
}

#define HASH_EDGE_DEF \
"void LG_HM_hash_edge (uint64_t *z, const uint64_t *x,\n"                      \
"    GrB_Index i, GrB_Index j, const uint64_t *seed)\n"                        \
"{\n"                                                                          \
"    uint64_t result = (i + j * i  + 0x9E3779B97F4A7C15LL) ;\n"                \
"    result = (result ^ (*x)) * 0xBF58476D1CE4E5B9LL ;\n"                      \
"    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9LL ;\n"            \
"    result = (result ^ (result >> 27)) * 0x94D049BB133111EBLL ;\n"            \
"    result = (result ^ (result >> 31)) ;\n"                                   \
"    (*z) = result ;\n"                                                        \
"}\n"

GrB_Info LAGraph_Hash_Matrix(
    uint64_t *hash,      // [output] hash
    const GrB_Matrix A,  // matrix to hash
    char *msg
) {
#if LG_SUITESPARSE_GRAPHBLAS_V10
    GrB_Matrix C = NULL;
    GrB_IndexUnaryOp lg_hash_edge = NULL;
    GrB_Index nrows, ncols;
    GRB_TRY (GrB_Matrix_nrows(&nrows, A)) ;
    GRB_TRY (GrB_Matrix_ncols(&ncols, A)) ;
    GRB_TRY (GrB_Matrix_new(&C, GrB_UINT64, nrows, ncols)) ;
    GRB_TRY (GxB_IndexUnaryOp_new(
        &lg_hash_edge, (GxB_index_unary_function) LG_HM_hash_edge,
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "LG_HM_hash_edge", HASH_EDGE_DEF)) ;
    
    // TODO: C takes extra memory which is not nessesary for this computation.
    // Compute without extra memory if possible.
    GRB_TRY (GrB_Matrix_apply_IndexOp_UINT64(
        C, NULL, NULL, lg_hash_edge, A, (uint64_t) 0, NULL));
    GRB_TRY (GrB_Matrix_reduce_UINT64(
        hash, GrB_BXOR_UINT64, GxB_BXOR_UINT64_MONOID, C, NULL)) ;
    LG_FREE_ALL;
    return GrB_SUCCESS;
#else
    return GrB_NOT_IMPLEMENTED;
#endif
}

GrB_Info LAGraph_Hash_Vector(
    uint64_t *hash,       // [output] hash
    const GrB_Vector v,   // Vector to hash
    char *msg
) {
#if LG_SUITESPARSE_GRAPHBLAS_V10
    return LAGraph_Hash_Matrix(hash, (GrB_Matrix) v, NULL);
#else
    return GrB_NOT_IMPLEMENTED;
#endif
}
