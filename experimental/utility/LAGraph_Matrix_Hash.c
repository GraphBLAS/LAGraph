//------------------------------------------------------------------------------
// LAGraph_Random_Matrix: generate a random matrix
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
    GrB_free(&x);                 \
    GrB_free(&s);                 \
    GrB_free(&lg_hash_edge);      \
    GrB_free(&lg_hash_edge_biop); \
    GrB_free(&lg_xor_hash);       \
}

static void hash_edge(
    uint64_t *z, 
    uint64_t *x,
    GrB_Index ix, 
    GrB_Index jx, 
    uint64_t *y,
    GrB_Index iy, 
    GrB_Index jy, 
    uint64_t thunk
) {
    uint64_t a = *x;
    a ^= ix + 0x9e3779b97f4a7c15ULL;
    a = (a << 13) | (a >> (64 - 13));
    a ^= jx + 0x9e3779b97f4a7c15ULL;
    a = (a << 17) | (a >> (64 - 17));
    *z = a;
}

GrB_Info LAGraph_Hash_Vector(uint64_t *hash, GrB_Vector v, char *msg){
    GrB_Vector x = NULL;
    GrB_Scalar s = NULL;
    GxB_IndexBinaryOp lg_hash_edge = NULL;
    GrB_BinaryOp lg_hash_edge_biop = NULL;
    GrB_Semiring lg_xor_hash = NULL;
    GrB_Index nrows;
    GRB_TRY (GrB_Vector_size(&nrows, v)) ;
    GRB_TRY (GrB_Scalar_new(&s, GrB_UINT64)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT64(s, 0)) ;
    GRB_TRY (GrB_Vector_new(&x, GrB_UINT64, nrows)) ;
    GRB_TRY (GrB_Vector_assign_UINT64(x, NULL, NULL, (uint64_t) 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GxB_IndexBinaryOp_new(&lg_hash_edge, (GxB_index_binary_function) hash_edge,
        GrB_UINT64, GrB_UINT64, GrB_UINT64, GrB_UINT64, NULL, NULL)) ;
    GRB_TRY (GxB_BinaryOp_new_IndexOp(&lg_hash_edge_biop, lg_hash_edge, s)) ;
    GRB_TRY (GrB_Semiring_new (&lg_xor_hash, GxB_BXOR_UINT64_MONOID, lg_hash_edge_biop)) ;
    GRB_TRY (GrB_mxv((GrB_Vector) s, NULL, NULL, lg_xor_hash, (GrB_Matrix) v, x, GrB_DESC_T0)) ;
    GRB_TRY (GrB_Scalar_extractElement_UINT64(hash, s)) ;
    LG_FREE_ALL;
    return GrB_SUCCESS;
}