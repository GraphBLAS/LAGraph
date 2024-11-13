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
    GrB_free (&P) ;                             \
    GrB_free (&M) ;                             \
    GrB_free (&M_fours) ;                       \
    GrB_free (&r_exists) ;                      \
    GrB_free(&new_hashed_edges);                \
    GrB_free(&buckets);                         \
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
    GrB_free (&r_pairs) ;                       \
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
    LAGraph_Free((void**)&indices, msg) ;       \
    LAGraph_Free((void**)&ramp, msg) ;          \
    LAGraph_Free((void**)&half_ramp, msg) ;     \
    LAGraph_Free((void**)&swap_type, msg) ;     \
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

    // swaps x e
    // Selected pairs for next batch of swaps
    // Each row contains 2 entries for the edges involved in a swap.
    GrB_Matrix P = NULL;

    // swaps x 4
    // Each row contains 4 entries corresponding to the verticies 
    // that are involved in the swap.
    GrB_Matrix M = NULL;
    GrB_Matrix M_fours = NULL; // M with exactly 4 entries

    GrB_Index n = 0, e = 0;

    // n x n 
    // Lower triangle of adjacency matrix
    GrB_Matrix A_tril = NULL ;

    // e x 1 random vectors
    GrB_Vector random_v = NULL, r_permute = NULL;

    //indicies for A
    GrB_Index *indices = NULL;

    // [0,2,0,3,0,3,0 . . .] numSwaps?
    GrB_Vector swapVals = NULL;

    // e x 2 Matrix that picks the edges for which we will swap values.
    GrB_Matrix swapMask = NULL;

    // Swaps 0,3 and 2,1
    GrB_Matrix swap_p = NULL;

    // This ramp will likely get recycled a few times. 
    GrB_Vector ramp_v = NULL;
    GrB_Vector hramp_v = NULL;

    // [0, ... , e]
    GrB_Index *ramp = NULL ;
    // [0,0,1,1,2,2,...]
    GrB_Index *half_ramp = NULL ;

    // Arrays to unpack edge permutation
    GrB_Index *edge_perm = NULL ;
    uint8_t *swap_type = NULL ;
    bool iso = false;

    // Reduced Vectors
    GrB_Vector r_exists = NULL, r_pairs;
    GrB_Matrix existMask = NULL;
    

    // Number of values kept in each phase
    GrB_Index n_keep;
    GrB_Index *arr_keep = NULL;
    void *junk = NULL;


    // swaps x 2 matrix which holds the hashes of each planned edge. 
    GrB_Matrix new_hashed_edges = NULL;

    // 2^60 x swaps holds the swaps in hash buckets
    GrB_Matrix buckets = NULL;

    GrB_Matrix p_buckets = NULL;
    
    // swaps used for "outdegree"
    GrB_Vector big_dense; 

    // e holds hashes of old edges
    GrB_Vector hashed_edges = NULL;

    // 2^60 holds the buckets in which hashes collided.
    GrB_Vector exists = NULL; 
    GrB_Index *hash_vals_new = NULL, *hash_vals = NULL;
    
    GrB_UnaryOp first_bit = NULL;

    // z = h_y(x)
    GrB_BinaryOp hash_seed = NULL;

    // [^],[h_y(x)]
    GrB_Semiring bxor_hash = NULL;

    GrB_Semiring bxor_first = NULL;

    GrB_Vector sort_h = NULL;
    GrB_Vector r_60;

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
    
    //----------------------------------------------------------- Init Operators
    GRB_TRY (GxB_UnaryOp_new (
        &first_bit, (GxB_unary_function) (&first_bit_equals),
        GrB_BOOL, GrB_UINT64, "first_bit_equals", FIRST_BIT_EQ
    )) ;

    GRB_TRY(GxB_BinaryOp_new(
        &hash_seed, (GxB_binary_function) (&hash_node),
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "hash_node", HASH_ONE
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
    
    // Extract adjacency matrix to make incidence matrix - E
    // get just the lower triangular entries
    //TODO: should I remove the diagonal? change the 0 if so
    // Arrays to extract A into

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
    GxB_Matrix_fprint(E, "E", GxB_SHORT, stdout);
    // Init Ramps --------------------------------------------------------------
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, e + 1)) ;
    GRB_TRY (GrB_Vector_new(&hramp_v, GrB_UINT64, e + 1)) ;
    GRB_TRY (GrB_Vector_new(&swapVals, GrB_BOOL, e)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (ramp_v, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
        hramp_v, NULL, NULL, GrB_DIV_UINT64, ramp_v, 2UL, NULL)) ;

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
    GRB_TRY (GrB_Vector_new (&r_pairs, GrB_INT64, e)) ;

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
        GRB_TRY (GrB_Matrix_new(&P, GrB_UINT8, e, e)) ; 
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

        

        GRB_TRY(GrB_Vector_apply(
            swapVals, NULL, NULL, first_bit, random_v, NULL
        )) ;
        

        GrB_Index perm_size;
        GRB_TRY (GxB_Vector_unpack_Full(
            r_permute, (void **)&edge_perm, &perm_size, &iso, NULL
        )) ;
        LG_ASSERT(!iso, GrB_NOT_IMPLEMENTED);

        
        // TODO: Should there be a function that takes in A matrix and computes
        // in or out degree 
        // GRB_TRY (GxB_Matrix_build_Scalar(
        //     P, ramp, edge_perm, one8, e
        // ));
        
        GRB_TRY (GxB_Matrix_pack_CSR(
            P, &ramp, &edge_perm, (void**) &val_of_P, (e + 1) * sizeof(GrB_Index),
            perm_size, sizeof(GrB_Index) , true, false, NULL
        ));

        // Pair edges. Take a random permutation and pair adjacent values.
        GRB_TRY (GrB_mxm(E, NULL, NULL, GxB_ANY_SECOND_UINT64, P, E, NULL)) ;
        
        GrB_Index arr_size, junk_size;
        GRB_TRY (GxB_Matrix_unpack_CSR(
            P, &ramp, &edge_perm, (void**) &val_of_P, &arr_size,
            &perm_size,  &junk_size, &iso, false, NULL
        ));

        //increase width of sorted so it can be used as a mask.
        GRB_TRY (GrB_mxm (swapMask, NULL, NULL, GxB_ANY_FIRST_BOOL,
            (GrB_Matrix) swapVals, (GrB_Matrix) dense_hash, GrB_DESC_T1)) ; 
        //swap vertexes in E randomly.
        GRB_TRY (GrB_mxm(
            E, swapMask, NULL, GxB_ANY_FIRST_UINT64, E, y, NULL)) ;
        
        GrB_Index E_bounds[3] = {swaps_per_loop * 2, e - swaps_per_loop * 2, 2};
        GRB_TRY (GrB_Matrix_new(E_split, GrB_UINT64, E_bounds[0], 2));
        GRB_TRY (GrB_Matrix_new(E_split + 1, GrB_UINT64, E_bounds[1], 2));
        
        GRB_TRY (GxB_Matrix_split(
            E_split, 2, 1, E_bounds, E_bounds + 2, E, NULL));

        M = E_split[0];
        GRB_TRY (GxB_Matrix_reshape(M, false, swaps_per_loop, 4, NULL)) ;

        GRB_TRY (GrB_Matrix_new(&new_hashed_edges, GrB_UINT64, swaps_per_loop, 2)) ;
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
        // GxB_Matrix_fprint(hash_s, "hash_s", GxB_SHORT, stdout);
        GRB_TRY(GrB_mxm(
            new_hashed_edges, NULL, NULL, bxor_hash, 
            M, hash_s, NULL
        ));
        GRB_TRY(GrB_mxv(
            hashed_edges, NULL, NULL, bxor_hash, 
            E, dense_hash, NULL
        ));

        // I will unpack and then reconstruct with hash as index.
        GRB_TRY(GxB_Matrix_unpack_FullR(
            new_hashed_edges, (void **) &hash_vals_new, &junk_size, &iso, NULL
        )) ;
        GRB_TRY(GxB_Vector_unpack_Full(
            hashed_edges, (void **) &hash_vals, &junk_size, &iso, NULL
        )) ;

        GRB_TRY (GrB_Matrix_new(
            &buckets, GrB_UINT64, 1ULL << 60, swaps_per_loop)) ;
        // Build hash buckets
        GRB_TRY(GxB_Vector_build_Scalar(
            exists, hash_vals, one64, e
        )) ;
        GRB_TRY(GxB_Matrix_build_Scalar(
            buckets, hash_vals_new, half_ramp, one64, swaps_per_loop * 2
        )) ;

        GRB_TRY (GrB_mxv(
            exists, NULL, GrB_PLUS_UINT64, LAGraph_plus_one_uint64, buckets, 
            big_dense, NULL
        )) ;

        GRB_TRY(GrB_Vector_select_UINT64(
            exists, NULL, NULL, GrB_VALUEGT_UINT64, exists, 1ull, GrB_DESC_R
        )) ;
        GRB_TRY(GrB_Vector_setElement_UINT64(exists, 1ll, (GrB_Index)0)) ;
        GRB_TRY(GrB_Vector_new(&r_exists, GrB_UINT8, swaps_per_loop)) ;
        
        // Search through array for bad swaps.
        GRB_TRY(GrB_vxm(
            r_exists, NULL, NULL, LAGraph_any_one_uint8, exists, buckets, NULL
        )) ;

        GRB_TRY (GrB_Vector_clear(exists)) ;

        // GxB_Vector_fprint(r_exists,"r_exists",GxB_SHORT, stdout);
        GRB_TRY (GrB_Vector_assign_INT8(
            r_exists, r_exists, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC)) ;
        GRB_TRY (GxB_Vector_unpack_CSC(
            r_exists, &arr_keep, &junk, &arr_size, &junk_size, &iso, &n_keep,
            NULL, NULL
        ));
        GRB_TRY(GrB_Matrix_new(&M_fours, GrB_UINT64, n_keep, 4)) ;

        GRB_TRY (GrB_Matrix_extract(
            M_fours, NULL, NULL, M, arr_keep, n_keep, 
            GrB_ALL, 0, NULL)) ;
        LG_TRY (LAGraph_Free(&junk, msg)) ;
        GRB_TRY (GrB_mxm(M_fours, NULL, NULL, GxB_ANY_FIRST_UINT64, M_fours, swap_p, NULL)) ;
        GRB_TRY (GrB_assign(M, NULL, NULL, M_fours, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GxB_Matrix_reshape(M, false, swaps_per_loop * 2, 2, NULL));
        
        GRB_TRY(GxB_Matrix_concat(E, E_split, 2, 1, NULL));

        // Free Matricies that have to be rebuilt
        FREE_LOOP ;

        num_attempts += swaps_per_loop;
        num_swaps += n_keep;
        swaps_per_loop = (n_keep * 3) / 2;
        swaps_per_loop = LAGRAPH_MAX(swaps_per_loop, 16) ;
        swaps_per_loop = LAGRAPH_MIN(swaps_per_loop, e / 3) ;

        // Maintain random Vector
        LG_TRY (LAGraph_Random_Next(random_v, msg)) ;
        printf("#####Made %ld swaps. Total %ld out of %ld. Attempting %ld swaps next.#####\n\n", n_keep, num_swaps, e * Q, swaps_per_loop);
    } 
    GrB_Index ind_size = 0;
    GRB_TRY (GxB_Matrix_unpack_FullC(
        E, (void **)&indices, &ind_size, &iso, NULL)) ;
    GRB_TRY (GxB_Matrix_build_Scalar(*A_new, indices, indices + e, one8, e));
    GRB_TRY (GrB_eWiseAdd(*A_new, NULL, NULL, GrB_PLUS_UINT8, *A_new,*A_new, GrB_DESC_T0)) ;
    GxB_Matrix_fprint(*A_new,"A_new",GxB_SHORT, stdout);
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
