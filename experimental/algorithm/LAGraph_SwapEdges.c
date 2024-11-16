//------------------------------------------------------------------------------
// LAGraph_SwapEdges: Randomly Swaps edges in a graph
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

#define FREE_LOOP                               \
{                                               \
    GrB_free (&M) ;                             \
    GrB_free (&M_fours) ;                       \
    GrB_free(&new_hashed_edges);                \
    GrB_free(&hashed_edges);                    \
    GrB_free (&big_dense) ;                     \
    LAGraph_Free((void **) &hash_vals_new, NULL);\
    LAGraph_Free((void **) &hash_vals, NULL);   \
}

#define LG_FREE_WORK                            \
{                                               \
    /* free any workspace used here */          \
    GrB_free (&E) ;                             \
    GrB_free (&A_tril) ;                        \
    GrB_free (&random_v) ;                      \
    GrB_free (&r_permute) ;                     \
    GrB_free (&ramp_v) ;                        \
    GrB_free (&hramp_v) ;                       \
    GrB_free (&swapVals) ;                      \
    GrB_free(&dense_hash);                      \
    GrB_free (&swap_p) ;                        \
    GrB_free (&first_bit) ;                     \
    GrB_free (&bxor_first) ;                    \
    GrB_free(&hash_s);                          \
    GrB_free (&hash_seed) ;                     \
    GrB_free (&r_60) ;                     \
    GrB_free(&exists);                          \
    GrB_free (&bxor_hash) ;                     \
    LAGraph_Free((void**)&dup_swaps, msg) ;       \
    LAGraph_Free((void**)&indices, msg) ;       \
    LAGraph_Free((void**)&ramp, msg) ;          \
    LAGraph_Free((void**)&half_ramp, msg) ;     \
    LAGraph_Free((void**)&edge_perm, msg) ;     \
    LAGraph_Free((void**) &val_of_P, msg);      \
    FREE_LOOP ;                                 \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free(A_new) ;                       \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

void first_bit_equals 
    (bool *z, const uint64_t *x)
    {
        (*z) = (bool) ((*x) & 1);
    }
#define FIRST_BIT_EQ                                                            \
"void first_bit_equals                                                       \n"\
"   (bool *z, const uint64_t *x)                                             \n"\
"   {                                                                        \n"\
"       (*z) = (bool) ((*x) & 1);                                            \n"\
"   }"

// creates [0,3,1,2,1,3,. . .] pattern from random vector.
void swap_pattern 
    (uint8_t *z, const uint64_t *x, int64_t i, int64_t j, const uint8_t *y)
    {
        (*z) = (uint8_t) (((i & 1) * 2) | (*x & 1));
    }
#define SWAP_PAT                                                                \
"void swap_pattern"                                                             \
    "(uint8_t *z, const uint64_t *x, int64_t i, int64_t j, const uint8_t *y)"   \
    "{"                                                                         \
        "(*z) = (uint8_t) (((i & 1) * 2) | (*x & 1));"                          \
    "}"

//Hashes any node with a simple Multiply-shift from 
// https://arxiv.org/pdf/1504.06804
// QUESTION: this hash is a bit simple but I doubt it will result in a ton of 
// collisions unless the input graph is very specifically constucted

void hash_node 
    (uint64_t *z, const uint64_t *x, const uint64_t *y)
{
    (*z) = ((*y) * (*x)) >> (4);
}
#define HASH_ONE                                                                \
"void hash_node                                                              \n"\
"    (uint64_t *z, const uint64_t *x, const uint64_t *y)                      \n"\
"{                                                                           \n"\
"    (*z) = ((*y) * (*x)) >> (4);                                            \n"\
"}"

void log_duplicate
    (uint64_t *z, const uint64_t *x, const uint64_t *y)
{
    **((int8_t **)y) = (int8_t) 0;
    *z = *x;
}
#define LOG_DUPLICATE                                                           \
"void log_duplicate                                                          \n"\
"    (uint64_t *z, const uint64_t *x, const uint64_t *y)                     \n"\
"{                                                                           \n"\
"    **((int8_t **)y) = (int8_t) 0;                                          \n"\
"    *z = *x;                                                                \n"\
"}"
int LAGraph_SwapEdges
(
    // output
    GrB_Matrix *A_new, //The adjacency matrix of G with edges randomly swapped
    // input: not modified
    LAGraph_Graph G,
    GrB_Index Q, // Swaps per edge
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL; // n x n Adjacency Matrix 

    // e x 2 with entries corresponding to verticies of an edge
    GrB_Matrix E = NULL, E_t = NULL; 

    // e entries. E_split[0] has those which are planning to swap.
    GrB_Matrix E_split[2] = {NULL, NULL}; 

    // swaps x 4
    // Each row contains 4 entries corresponding to the verticies 
    // that are involved in the swap.
    GrB_Matrix M = NULL;
    GrB_Matrix M_fours = NULL; // M with exactly 4 entries

    // n = |V| e = |E|
    GrB_Index n = 0, e = 0;

    // n x n 
    // Lower triangle of adjacency matrix
    GrB_Matrix A_tril = NULL ;

    // e x 1 random vectors
    GrB_Vector random_v = NULL, r_permute = NULL;

    // indicies for A
    GrB_Index *indices = NULL;

    // IDK
    GrB_Vector swapVals = NULL;

    // e x 2 Matrix that picks the edges for which we will swap values.
    GrB_Matrix swapMask = NULL;

    // Swaps 2,1
    GrB_Matrix swap_p = NULL;

    GrB_Vector ramp_v = NULL;
    GrB_Vector hramp_v = NULL;

    // [0, ... , e]
    GrB_Index *ramp = NULL ;
    // [0,0,1,1,2,2,...]
    GrB_Index *half_ramp = NULL ;

    // Arrays to unpack edge permutation
    GrB_Index *edge_perm = NULL ;
    bool iso = false;

    // Number of values kept in each phase
    GrB_Index n_keep;
    GrB_Index *arr_keep = NULL;
    void *junk = NULL;


    // swaps x 2 matrix which holds the hashes of each planned edge. 
    GrB_Matrix new_hashed_edges = NULL;

    // GrB_Matrix p_buckets = NULL;
    
    // swaps used for "outdegree"
    GrB_Vector big_dense; 

    // e holds hashes of old edges
    GrB_Vector hashed_edges = NULL;

    // 2^60 holds the buckets in which hashes collided.
    GrB_Vector exists = NULL; 
    GrB_Vector new_edges_h = NULL; 
    GrB_Index *hash_vals_new = NULL, *hash_vals = NULL;
    
    GrB_UnaryOp first_bit = NULL;

    // z = h_y(x)
    GrB_BinaryOp hash_seed = NULL;

    // [^],[h_y(x)]
    GrB_Semiring bxor_hash = NULL;

    GrB_Semiring bxor_first = NULL;

    GrB_BinaryOp duplicate = NULL;

    int8_t *dup_swaps = NULL;
    GrB_Vector dup_swaps_v = NULL;
    GrB_Vector not_pointers = NULL;
    uint64_t *not_ptrs = NULL; 

    GrB_Vector sort_h = NULL;
    GrB_Vector r_60 = NULL;

    // Constants ---------------------------------------------------------------

    uint64_t *val_of_P = NULL;

    // 4 x 2 Matrix Hashes cols 0,2 and 1,3 
    GrB_Matrix hash_s = NULL; 

    // 2 x 2 Used to shuffle edges [[0 1],[1 0]]
    GrB_Matrix y = NULL;

    // 2 x 1 with hash seed as value.
    GrB_Vector dense_hash = NULL;

    GrB_Scalar zero8 = NULL, one8 = NULL, one64 = NULL ;
    
    //--------------------------------------------------------------------------
    // Check inputs TODO
    //--------------------------------------------------------------------------
    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED,
        LAGRAPH_INVALID_GRAPH, 
        "G must be undirected"
    ) ;

    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, 
        "G->nself_edges must be zero") ;
    // GRB_TRY (GrB_get(A, (void*)&type, GrB_EL_TYPE_CODE)) ;
    // LG_ASSERT_MSG (type == GrB_BOOL_CODE, LAGRAPH_INVALID_GRAPH, 
    //     "A must be type boolean") ;

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;    
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(A_new, GrB_UINT8, n, n)) ;

    // Extract lower triangular edges.
    GRB_TRY (GrB_Matrix_new (&A_tril, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_nvals(&e, A_tril)) ;
    GRB_TRY (GrB_Matrix_new(&E, GrB_UINT64, e, 2)) ;
    GRB_TRY (GrB_Matrix_new(&E_t, GrB_UINT64, 2, e)) ;
    
    //Init Operators -----------------------------------------------------------
    GRB_TRY (GxB_UnaryOp_new (
        &first_bit, (GxB_unary_function) (&first_bit_equals),
        GrB_BOOL, GrB_UINT64, "first_bit_equals", FIRST_BIT_EQ
    )) ;

    GRB_TRY(GxB_BinaryOp_new(
        &hash_seed, (GxB_binary_function) (&hash_node),
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "hash_node", HASH_ONE
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &duplicate, (GxB_binary_function) (&log_duplicate),
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "log_duplicate", LOG_DUPLICATE
    )) ;
    // I use a bit wise xor to combine the hashes since the same column number 
    // will not appear twice in my multiplication and I want combination to be 
    // commutative.
    GRB_TRY(GrB_Semiring_new(
        &bxor_hash, GxB_BXOR_UINT64_MONOID, hash_seed 
    )) ;
    GRB_TRY(GrB_Semiring_new(
        &bxor_first, GxB_BXOR_UINT64_MONOID , GrB_FIRST_UINT64 
    )) ;
    // count swaps 
    GrB_Index num_swaps = 0, num_attempts = 0, swaps_per_loop = e / 3 ;

    // Make E Matrix -----------------------------------------------------------
    LG_TRY (LAGraph_Malloc ((void**)(&indices), 2ull * e, sizeof(GrB_Index), msg)) ;
    GRB_TRY (
        GrB_Matrix_extractTuples_BOOL (indices, indices + e, NULL, &e, A_tril)
        ) ;
    GRB_TRY (GxB_Matrix_pack_FullR (
        E_t, (void **)&indices, 2ull * e * sizeof(GrB_Index), false, NULL
    )) ;
    GRB_TRY (GrB_transpose(E, NULL, NULL, E_t, NULL));
    GrB_free(&E_t);
    GRB_TRY (GrB_Vector_new(&exists, GrB_UINT64, 1ULL << 60)) ;
    // GxB_Matrix_fprint(E, "E", GxB_SHORT, stdout);
    // Init Ramps --------------------------------------------------------------
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, e + 1)) ;
    GRB_TRY (GrB_Vector_new(&hramp_v, GrB_UINT64, e + 1)) ;
    GRB_TRY (GrB_Vector_new(&swapVals, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (ramp_v, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
        hramp_v, NULL, NULL, GrB_DIV_UINT64, ramp_v, 2UL, NULL)) ;
    GRB_TRY(GrB_Vector_dup(&not_pointers, hramp_v)) ;
    GrB_Index ramp_size;
    GRB_TRY (GxB_Vector_unpack_Full (
        ramp_v, (void **)&ramp, &ramp_size, &iso, NULL)) ;
    GRB_TRY (GxB_Vector_unpack_Full (
        hramp_v, (void **)&half_ramp, &ramp_size, &iso, NULL)) ;

    // Init Constants ----------------------------------------------------------
    GRB_TRY (GrB_Scalar_new (&zero8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_new (&one8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_new (&one64, GrB_UINT64)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (zero8, 0)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (one8, 1)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT64 (one64, 1ull)) ;

    GRB_TRY (GrB_Vector_new (&dense_hash, GrB_UINT64, 2)) ;
    GRB_TRY (GrB_Matrix_new (&y, GrB_UINT8, 2, 2)) ;
    GRB_TRY (GrB_Matrix_new (&hash_s, GrB_UINT64, 4, 2)) ;
    GRB_TRY (GrB_Matrix_new (&swap_p, GrB_UINT8, 4, 4)) ;
    GRB_TRY (GrB_Vector_new (&new_edges_h, GrB_UINT64, 1ull << 60)) ;

    GRB_TRY (GrB_Vector_assign_UINT8 (
        dense_hash, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_assign_UINT8 (
        y, NULL, NULL, 0, GrB_ALL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY(GrB_Matrix_select_UINT64(y, NULL, NULL, GrB_OFFDIAG, y, 0, NULL)) ;

    LG_TRY (LAGraph_Malloc ((void**)(&val_of_P), 1, sizeof(GrB_Index), msg)) ;
    val_of_P[0] = 1;

    GrB_Index hcols[] = {0, 0, 1, 1};
    GrB_Index srows[] = {0, 1, 2, 3};
    GrB_Index scols[] = {0, 2, 1, 3};
    GRB_TRY(GxB_Matrix_build_Scalar(
        hash_s, scols, hcols, one64, 4
    )) ;
    GRB_TRY(GxB_Matrix_build_Scalar(
        swap_p, srows, scols, one64, 4
    )) ;

    // Make Random -------------------------------------------------------------
    GRB_TRY (GrB_Vector_new(&random_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_60, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        random_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;
    //TODO: Change seed
    LG_TRY(
        LAGraph_Random_Seed(random_v, 1548945616ul, msg));

    printf("Entering loop, Good Luck:\n") ;
    while(num_swaps < e * Q && num_attempts < e * Q * 5)
    {
        // Coming into the loop: 
        // E must be the incidence matrix of the new graph. W/o self edges nor 
        // parallel edges. Each row must have exactly two distinct values.
        // random_v has a radom dense vector.
        GRB_TRY (GrB_Matrix_new (&swapMask, GrB_BOOL, e, 2)) ;
        // GRB_TRY (GrB_Matrix_new (&p_buckets, GrB_UINT8, e, 2)) ;


        GRB_TRY (GxB_Vector_sort (
            NULL, r_permute, GrB_LT_UINT64, random_v, GrB_NULL
        )) ;
        // GrB_Index *rand_arr, rand_size;
        // GRB_TRY(GrB_Vector_apply_BinaryOp1st_UINT64(
        //     r_60, NULL, NULL, GrB_BAND_UINT64, (uint64_t) 0xFFFFFFFFFFFFFFF, random_v, NULL
        // )) ;
        // GRB_TRY (GxB_Vector_unpack_Full(
        //     r_60, (void **)&rand_arr, &rand_size, &iso, NULL
        // )) ;
        // Hash (rand_arr)
        // Matrix_Build (rand_arr, ramp, ramp, )

        

        GRB_TRY(GrB_Vector_apply(
            swapVals, NULL, NULL, first_bit, random_v, NULL
        )) ;
        

        GrB_Index perm_size, arr_size, junk_size;
        GRB_TRY (GxB_Vector_unpack_Full(
            r_permute, (void **)&edge_perm, &perm_size, &iso, NULL
        )) ;
        LG_ASSERT(!iso, GrB_NOT_IMPLEMENTED);

        GRB_TRY (GrB_Matrix_extract(
            E, NULL, NULL, E, edge_perm, e, GrB_ALL, 0, NULL
        )) ;
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            swapVals, swapVals, NULL, GxB_BXOR_UINT64_MONOID, E, NULL));
        GrB_Matrix xor_diag = NULL;
        GRB_TRY (GrB_Matrix_diag(&xor_diag, swapVals, 0));
        GRB_TRY (GrB_mxm(
            E, NULL, NULL, GxB_BXOR_BXOR_UINT64, xor_diag, E, NULL)) ;
        //increase width of sorted so it can be used as a mask.
        // GRB_TRY (GrB_mxm (swapMask, NULL, NULL, GxB_ANY_FIRST_BOOL,
        //     (GrB_Matrix) swapVals, (GrB_Matrix) dense_hash, GrB_DESC_T1)) ; 
        //swap vertexes in E randomly.
        // GRB_TRY (GrB_mxm(
        //     E, swapMask, NULL, GxB_ANY_FIRST_UINT64, E, y, NULL)) ;
        
        GrB_Index E_bounds[3] = {swaps_per_loop * 2, e - swaps_per_loop * 2, 2};
        GRB_TRY (GrB_Matrix_new(E_split, GrB_UINT64, E_bounds[0], 2));
        GRB_TRY (GrB_Matrix_new(E_split + 1, GrB_UINT64, E_bounds[1], 2));
        
        GRB_TRY (GxB_Matrix_split(
            E_split, 2, 1, E_bounds, E_bounds + 2, E, NULL));

        M = E_split[0];
        GRB_TRY (GxB_Matrix_reshape(M, false, swaps_per_loop, 4, NULL)) ;

        // Hash Edges ----------------------------------------------------------
        GRB_TRY (GrB_Matrix_new(
            &new_hashed_edges, GrB_UINT64, swaps_per_loop, 2)) ;
        GRB_TRY (GrB_Vector_new(&hashed_edges, GrB_UINT64, e)) ;
        GRB_TRY (GrB_Vector_new(&big_dense, GrB_UINT64, swaps_per_loop)) ;

        //This Scalar is used in the hash funtion and should be something near 
        // 2^60 and odd.
        GRB_TRY (GrB_Matrix_assign_UINT64(
            hash_s, hash_s, NULL, 0xFB21C651E98DF25ULL, 
            GrB_ALL, 0, GrB_ALL, 0, GrB_DESC_S
        )) ;
        GRB_TRY (GrB_Vector_assign_UINT64(
            dense_hash, NULL, NULL, 0xFB21C651E98DF25ULL, GrB_ALL, 0, NULL
        )) ;
        GRB_TRY (GrB_Vector_assign_UINT64(
            big_dense, NULL, NULL, (uint64_t) 0, GrB_ALL, 0, NULL
        )) ;
        GRB_TRY(GrB_mxm(
            new_hashed_edges, NULL, NULL, bxor_hash, 
            M, hash_s, NULL
        )) ;
        GRB_TRY(GrB_mxv(
            hashed_edges, NULL, NULL, bxor_hash, 
            E, dense_hash, NULL
        ));

        // I will unpack and then reconstruct with hash as index.
        GRB_TRY(GxB_Matrix_unpack_BitmapR(
            new_hashed_edges, &dup_swaps, (void **) &hash_vals_new, 
            &junk_size, &junk_size, &iso, &junk_size, NULL
        )) ;
        GRB_TRY(GxB_Vector_unpack_Full(
            hashed_edges, (void **) &hash_vals, &junk_size, &iso, NULL
        )) ;
        GRB_TRY(GrB_Vector_apply_BinaryOp1st_UINT64(
            not_pointers, NULL,NULL, GrB_PLUS_UINT64, (uint64_t) dup_swaps, 
            not_pointers, NULL
        )) ;
        GRB_TRY(GxB_Vector_unpack_Full(
            not_pointers, (void **) &not_ptrs, &arr_size, &iso, NULL
        )) ;

        // Build Hash Buckets --------------------------------------------------
        GRB_TRY(GrB_Vector_build_UINT64(
            new_edges_h, hash_vals_new, not_ptrs, swaps_per_loop * 2, duplicate
        )) ;
        GRB_TRY(GxB_Vector_build_Scalar(
            exists, hash_vals, one64, e
        )) ;
        GRB_TRY(GrB_Vector_setElement_UINT64(exists, 1ull, (GrB_Index)0)) ;
        GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(
            new_edges_h, NULL, NULL, duplicate, exists, new_edges_h, NULL
        ));

        // TODO: not sure I need this
        GRB_TRY(GrB_wait(new_edges_h, GrB_COMPLETE)) ;
        // GRB_TRY(GrB_Vector_new(
        //     &dup_swaps_v, GrB_BOOL, swaps_per_loop * 2
        // )) ;
        GRB_TRY(GxB_Vector_pack_Full(
            not_pointers, (void **) &not_ptrs, arr_size, iso, NULL
        )) ;
        GRB_TRY(GrB_Vector_apply_BinaryOp2nd_UINT64(
            not_pointers, NULL,NULL, GrB_MINUS_UINT64, not_pointers,
            (uint64_t) dup_swaps, NULL
        )) ;
        LG_TRY(LAGraph_Malloc(
            (void **) &arr_keep, swaps_per_loop, sizeof(uint64_t), msg) ;)
        n_keep = 0;
        for (int64_t i = 0; i < swaps_per_loop; i++)
            if(dup_swaps[i])
                arr_keep[n_keep++] = i;
        GRB_TRY (GrB_Vector_clear(exists)) ;
        GRB_TRY (GrB_Vector_clear(new_edges_h)) ;
        // bool *iso_dup_value = NULL;
        // LAGraph_Malloc((void **) &iso_dup_value, 1, 1, msg) ;
        // iso_dup_value[0] = true;
        // GRB_TRY(GxB_Vector_pack_Bitmap(
        //     dup_swaps_v, &dup_swaps, (void **)&iso_dup_value, 
        //     swaps_per_loop * 2, 1, true, swaps_per_loop * 2, NULL
        // ));
        // GxB_Vector_fprint(dup_swaps_v, "Dups", GxB_SHORT, stdout) ;
        // Search through array for bad swaps.

        // Swap Good Edges -----------------------------------------------------
        GRB_TRY(GrB_Matrix_new(&M_fours, GrB_UINT64, n_keep, 4)) ;

        GRB_TRY (GrB_Matrix_extract(
            M_fours, NULL, NULL, M, arr_keep, n_keep, GrB_ALL, 0, NULL
        )) ;
        GRB_TRY (GrB_mxm(
            M_fours, NULL, NULL, GxB_ANY_FIRST_UINT64, M_fours, swap_p, NULL
        )) ;
        GRB_TRY (GrB_assign(
            M, NULL, NULL, M_fours, arr_keep, n_keep, GrB_ALL, 0, NULL
        )) ;
        GRB_TRY (GxB_Matrix_reshape(M, false, swaps_per_loop * 2, 2, NULL));
        GRB_TRY(GxB_Matrix_concat(E, E_split, 2, 1, NULL));

        
        FREE_LOOP ; // Free Matricies that have to be rebuilt

        // Adjust number of swaps to do next.
        num_attempts += swaps_per_loop;
        num_swaps += n_keep;
        swaps_per_loop = (n_keep * 3) / 2;
        swaps_per_loop = LAGRAPH_MAX(swaps_per_loop, 16) ;
        swaps_per_loop = LAGRAPH_MIN(swaps_per_loop, e / 3) ;

        LG_TRY (LAGraph_Random_Next(random_v, msg)) ;
        printf("#####Made %ld swaps. Total %ld out of %ld. Attempting %ld swaps next.#####\n\n", n_keep, num_swaps, e * Q, swaps_per_loop);
    } 
    // Build Output Matrix
    GrB_Index ind_size = 0;
    GRB_TRY (GxB_Matrix_unpack_FullC(
        E, (void **)&indices, &ind_size, &iso, NULL)) ;
    GRB_TRY (GxB_Matrix_build_Scalar(*A_new, indices, indices + e, one8, e));
    GRB_TRY (GrB_eWiseAdd(
        *A_new, NULL, NULL, GrB_PLUS_UINT8, *A_new,*A_new, GrB_DESC_T0
    )) ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
