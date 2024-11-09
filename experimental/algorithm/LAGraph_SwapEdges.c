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
    GrB_free (&P) ;                         \
    GrB_free (&M) ;                             \
    GrB_free (&M_fours) ;                       \
    GrB_free (&pairs_4s) ;                      \
    GrB_free (E_arranged) ;                     \
    GrB_free (E_arranged + 1) ;                 \
    GrB_free (E_arranged + 2) ;                 \
    GrB_free (&pairs_new) ;                     \
    GrB_free (&M_outdeg) ;                      \
    GrB_free (&r_exists) ;                      \
    GrB_free(&new_hashed_edges);                \
    GrB_free(&buckets);                         \
    GrB_free(&hashed_edges);                    \
    GrB_free(&exists);                          \
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
    GrB_free (&r_sorted) ;                      \
    GrB_free (&x) ;                             \
    GrB_free (&r_pairs) ;                       \
    GrB_free (&ramp_v) ;                        \
    GrB_free (&hramp_v) ;                       \
    GrB_free (&swapVals) ;                      \
    GrB_free(&dense_hash);                     \
    GrB_free (&swap_p) ;                      \
    GrB_free (&first_bit) ;                     \
    GrB_free (&swap_op) ;                       \
    GrB_free (&bxor_first) ;                      \
    GrB_free(&hash_s);                      \
    GrB_free (&hash_one) ;                      \
    GrB_free (&hash_two) ;                      \
    GrB_free (&hash_seed) ;                     \
    GrB_free (&bhash_two) ;                     \
    GrB_free (&hash_edges) ;                    \
    GrB_free (&hash_edges2) ;                   \
    LAGraph_Free((void**)&row_indices, msg) ;   \
    LAGraph_Free((void**)&col_indices, msg) ;   \
    LAGraph_Free((void**)&ramp, msg) ;          \
    LAGraph_Free((void**)&half_ramp, msg) ;     \
    LAGraph_Free((void**)&swap_type, msg) ;     \
    LAGraph_Free((void**)&edge_perm, msg) ;     \
    LAGraph_Free((void**)&swaps, msg) ;         \
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

//hash that places the values for the pair of nodes corresponding to one edge in 
//one part of an array
void hash_pair 
    (uint64_t *z, const uint8_t *x, GrB_Index i_x, GrB_Index j_x, 
    const uint64_t *y, GrB_Index i_y, GrB_Index j_y, const int *nbits)
{
    (*z) = (((*x) & 1) ^ j_y) * ((*y) * j_x) >> (64 - *nbits);
}
#define HASH_TWO                                                                \
"void hash_pair                                                              \n"\
"    (uint64_t *z, const uint8_t *x, GrB_Index i_x, GrB_Index j_x,           \n"\
"    const uint64_t *y, GrB_Index i_y, GrB_Index j_y, const int *nbits)      \n"\
"{                                                                           \n"\
"    (*z) = (((*x) & 1) ^ j_y) * ((*y) * j_x) >> (64 - *nbits);              \n"\
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
    GrB_Matrix E = NULL; 

    // cmatrix [E_selected, M_1, M_2]
    GrB_Matrix E_arranged[3] = {NULL, NULL, NULL};

    // swaps x e
    // Selected pairs for next batch of swaps
    // Each row contains 2 entries for the edges involved in a swap.
    GrB_Matrix P = NULL, pairs_4s = NULL;
    GrB_Matrix pairs_new = NULL;

    // swaps x 4
    // Each row contains 4 entries corresponding to the verticies 
    // that are involved in the swap.
    GrB_Matrix M = NULL;
    GrB_Matrix M_xor = NULL; // swaps x 2 with a zero in self edges.
    GrB_Matrix M_fours = NULL; // M with exactly 4 entries

    GrB_Index n = 0, e = 0;

    // Copied vars from LAGraph_Incidence_Matrix
    GrB_Matrix A_tril = NULL ;

    // random vector 
    GrB_Vector random_v = NULL, r_permute = NULL, r_sorted = NULL;

    //indicies for A
    GrB_Index *row_indices = NULL, *col_indices = NULL ;

    // [0,2,0,3,0,3,0 . . .] numSwaps?
    GrB_Vector swapVals = NULL;

    // e x 2 Matrix that picks the edges for which we will swap values.
    GrB_Matrix swapMask = NULL;

    
    GrB_Matrix swap_p = NULL;

    // This ramp will likely get recycled a few times. 
    GrB_Vector ramp_v = NULL;
    GrB_Vector hramp_v = NULL;

    // C array [0, ... , e-1]
    GrB_Index *ramp = NULL ;
    GrB_Index *half_ramp = NULL ;

    // Arrays to unpack edge permutation
    GrB_Index *edge_perm = NULL ;
    uint8_t *swap_type = NULL ;
    bool iso = false;

    // Reduced Vectors
    GrB_Vector M_outdeg = NULL, r_exists = NULL, r_pairs;
    GrB_Matrix existMask = NULL;
    
    // 4 x 1 full vector 
    GrB_Vector x = NULL;

    // Number of values kept in each phase
    GrB_Index n_keep;
    GrB_Index *arr_keep = NULL;
    void *junk = NULL;

    // 4 x 2 Matrix Hashes cols 0,2 and 1,3 
    GrB_Matrix hash_s = NULL; 
    GrB_Matrix new_hashed_edges = NULL, buckets = NULL;
    GrB_Vector dense_hash = NULL, hashed_edges = NULL, exists = NULL; 
    GrB_Index *hash_vals_new = NULL, *hash_vals = NULL;
    uint8_t *swaps = NULL;


    GrB_UnaryOp first_bit = NULL;
    GrB_IndexUnaryOp swap_op = NULL;
    GzB_IndexBinaryOp hash_one = NULL, hash_two = NULL;
    GrB_BinaryOp bhash_two = NULL, hash_seed = NULL;
    GrB_Semiring hash_edges = NULL, hash_edges2 = NULL;

    GrB_Semiring bxor_first = NULL;

    A = G->A ;
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
    LAGRAPH_TRY(LAGraph_Random_Init(msg)) ;
    
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(A_new, GrB_UINT8, n, n)) ;

    // Extract lower triangular edges.
    GRB_TRY (GrB_Matrix_new (&A_tril, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_nvals(&e, A_tril)) ;

    GRB_TRY (GrB_Matrix_new(&E, GrB_UINT64, e, 2)) ;
    

    GRB_TRY (GxB_UnaryOp_new (
        &first_bit, (GxB_unary_function) (&first_bit_equals),
        GrB_BOOL, GrB_UINT64, "first_bit_equals", FIRST_BIT_EQ
    )) ;
    GRB_TRY (GxB_IndexUnaryOp_new (
       &swap_op, (GxB_index_unary_function) (&swap_pattern), 
       GrB_UINT8, GrB_UINT64, GrB_UINT8, "swap_pattern", SWAP_PAT
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &hash_seed, (GxB_binary_function) (&hash_node),
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "hash_node", HASH_ONE
    )) ;

    GRB_TRY(GzB_IndexBinaryOp_new(
        &hash_two, (GzB_index_binary_function) (&hash_pair), 
        GrB_UINT64, GrB_UINT8, GrB_UINT64, GrB_INT32, "hash_pair", HASH_TWO
    )) ;

    GrB_Scalar sizeTheta ; // hash size
    GrB_Scalar_new (&sizeTheta, GrB_INT32) ;
    GrB_Scalar_setElement_INT32 (sizeTheta, 60);
    GRB_TRY(GzB_BinaryOp_new_IndexOp(&bhash_two, hash_two, sizeTheta));
    
    // I use a bit wise xor to combine the hashes since the same column number 
    // will not appear twice in my multiplication and I want combination to be 
    // commutative.
    GRB_TRY(GrB_Semiring_new(
        &hash_edges, GxB_BXOR_UINT64_MONOID, hash_seed 
    )) ;
    GRB_TRY(GrB_Semiring_new(
        &hash_edges2, GxB_BXOR_UINT64_MONOID, bhash_two
    )) ;
    GRB_TRY(GrB_Semiring_new(
        &bxor_first, GxB_BXOR_UINT64_MONOID , GrB_FIRST_UINT64 
    )) ;

    GrB_Index num_swaps = 0, num_attempts = 0; 
    // Q: Should this decrease if edges are interfering alot?
    // Bound number of swaps by E-2 * (max deg)?
    // Yes might be good to decrease this if you get alot of intrf
    // or increase if intrf low
    // maybe a * #swaps that worked last iteration
    GrB_Index swaps_per_loop = e / 3 ; // Make this a cap
    
    

    // This is length e to let every edge have a chance to swap.
    // QUESTION: Perhaps it's better to just do random_v % e and deal with 
    // duplicate edges and self paired edges later. 
    // (inteference matrix will detect)

    GrB_Matrix y;
    GRB_TRY (GrB_Vector_new(&random_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_sorted, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new (&dense_hash, GrB_UINT64, 2)) ;
    GRB_TRY (GrB_Matrix_new (&y, GrB_UINT8, 2, 2)) ;
    GRB_TRY (GrB_Matrix_new (&hash_s, GrB_UINT64, 4, 2)) ;
    GRB_TRY (GrB_Matrix_new (&swap_p, GrB_UINT8, 4, 4)) ;
    GRB_TRY (GrB_Vector_new (&r_pairs, GrB_INT64, e)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_INT8, 4)) ;



    // Extract adjacency matrix to make incidence matrix - E
    // get just the lower triangular entries
    //TODO: should I remove the diagonal? change the 0 if so
    // Arrays to extract A into

    LG_TRY (LAGraph_Malloc ((void**)(&row_indices), 2 * e, sizeof(GrB_Index), msg)) ;
    GRB_TRY (
        GrB_Matrix_extractTuples_BOOL (row_indices, row_indices + e, NULL, &e, A_tril)
        ) ;
    GRB_TRY (GxB_Matrix_pack_FullC (
        E, (void **)&row_indices, 2 * e * sizeof(GrB_Index), false, NULL
    )) ;

    //Build large ramp
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&hramp_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&swapVals, GrB_BOOL, e)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (ramp_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    // [0,0,1,1,2,2,...]
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
        hramp_v, NULL, NULL, GrB_DIV_UINT64, ramp_v, 2UL, NULL)) ;

    GrB_Index ramp_size;
    GRB_TRY (GxB_Vector_unpack_Full (
        ramp_v, (void **)&ramp, &ramp_size, &iso, NULL)) ;
    GRB_TRY (GxB_Vector_unpack_Full (
        hramp_v, (void **)&half_ramp, &ramp_size, &iso, NULL)) ;

    GrB_Scalar zero8 = NULL, one8 = NULL, one64 = NULL ;
    GRB_TRY (GrB_Scalar_new (&zero8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_new (&one8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_new (&one64, GrB_UINT64)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (zero8, 0)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (one8, 1)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT64 (one64, 1ull)) ;


    GRB_TRY (GrB_Vector_assign_UINT8 (
        dense_hash, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_assign_UINT8 (
        y, NULL, NULL, 0, GrB_ALL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign_UINT8 (
        x, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY(GrB_Matrix_select_UINT64(y, NULL, NULL, GrB_OFFDIAG, y, 0, NULL)) ;

    GrB_Index hcols[] = {0, 0, 1, 1};
    GrB_Index srows[] = {0, 1, 2, 3};
    GrB_Index scols[] = {0, 2, 1, 3};
    GRB_TRY(GxB_Matrix_build_Scalar(
        hash_s, scols, hcols, one64, 4
    )) ;
    GRB_TRY(GxB_Matrix_build_Scalar(
        swap_p, srows, scols, one64, 4
    )) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        random_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;

    //TODO: Change seed
    LG_TRY(
        LAGraph_Random_Seed(random_v, 1548945616ul, msg));

    printf("Entering loop, Good Luck:\n") ;
    while(num_swaps < e * Q && num_attempts < e * Q * 5)
    {
        // GRB_TRY (GrB_mxm(
        //     *A_new, NULL, NULL, LAGraph_plus_one_uint8, E, E, GrB_DESC_RT0)) ;
        // GRB_TRY(GrB_Matrix_select_UINT8(
        //     *A_new, NULL, NULL, GrB_OFFDIAG, *A_new, (uint8_t) 0, NULL)) ;
        // GRB_TRY(GrB_Matrix_select_UINT8(
        //     *A_new, NULL, NULL, GrB_VALUEEQ_UINT8, *A_new, (uint8_t) 2, NULL)) ;
        // GRB_TRY(GrB_Matrix_nvals(&nvals, *A_new));
        // if(nvals)
        // {
        //     printf("SELF EDGE MADE \n\n");
        //     GxB_Matrix_fprint(*A_new,"Hashes",GxB_COMPLETE, stdout);
        //     GRB_TRY(-1);
        // }
        // Coming into the loop: 
        // E must be the incidence matrix of the new graph. W/o self edges nor 
        // parallel edges. Each row must have exactly one 0 and one 1. 
        // They do not have to be randomly assigned.
        // random_v has a radom dense vector.

        GRB_TRY (GrB_Matrix_new(&P, GrB_UINT8, e, e)) ; 
        GRB_TRY (GrB_Matrix_new (&M_xor, GrB_UINT8, swaps_per_loop, 2)) ;
        GRB_TRY (GrB_Matrix_new (&swapMask, GrB_BOOL, e, 2)) ;
        GRB_TRY (GrB_Vector_new (&M_outdeg, GrB_UINT8, swaps_per_loop)) ;
        
        
        GRB_TRY (GxB_Vector_sort (
            r_sorted, r_permute, GrB_LT_UINT64, random_v, GrB_NULL
        )) ;
        GRB_TRY(GrB_Vector_apply(swapVals, NULL, NULL, first_bit, 
            r_sorted, NULL)) ;

        // TODO: try and make this more efficient maybe just rand % e rather 
        // than a permutation
        GrB_Index perm_size;
        GRB_TRY (GxB_Vector_unpack_Full(
            r_permute, (void **)&edge_perm, &perm_size, &iso, NULL
        )) ;
        LG_ASSERT(!iso, GrB_NOT_IMPLEMENTED);

        // Pair edges. Take a random permutation and pair adjacent values.
        // pairs wants:
        // cols: [0, 0, 1, 1, 2, 2, 3, 3, . . .] (ramp / 2)
        // rows: [ random permutation 1 - nvals] (LAGraph_Random_Seed + sort)
        // vals: [1, 3, 0, 2, 0, 2, 0, 3, . . .]
            // Start with [0, 2, 0, 2, 0, 2, . . .]
            // for i from 1 to nvals by 2s:
            //      vals[i] += vals[i] > vals[i-1]
        // I know it will be alot less safe but is it more efficient to pack
        // pairs as needed? Doubt it.

        // uint64_t *val_of_P = NULL;
        // LG_TRY (LAGraph_Malloc ((void**)(&val_of_P), 1, sizeof(GrB_Index), msg)) ;
        // val_of_P[0] = 1;
        // GRB_TRY (GxB_Matrix_pack_CSR(
        //     P, &ramp, &edge_perm, (void**) &val_of_P, (e + 1) * sizeof(GrB_Index),
        //     perm_size, sizeof(GrB_Index) , true, false, NULL
        // ));
        // TODO: Should there be a function that takes in A matrix and computes
        // in or out degree 
        GRB_TRY (GxB_Matrix_build_Scalar(
            P, ramp, edge_perm, one8, e
        ));

        //Permute the edges
        GRB_TRY (GrB_mxm(E, NULL, NULL, GxB_ANY_SECOND_UINT64, P, E, NULL)) ;


        //increase width of sorted so it can be used as a mask.
        GRB_TRY (GrB_mxm (swapMask, NULL, NULL, GxB_ANY_FIRST_BOOL,
            (GrB_Matrix) swapVals, (GrB_Matrix) dense_hash, GrB_DESC_T1)) ; 
        //swap vertexes in E randomly.
        GRB_TRY (GrB_mxm(
            E, swapMask, NULL, GxB_ANY_FIRST_UINT64, E, y, NULL)) ;
        GRB_TRY (GrB_Matrix_dup(&M,E));
        GRB_TRY (GrB_Matrix_resize(M, swaps_per_loop * 2, 2));
        GRB_TRY (GxB_Matrix_reshape(M, false, swaps_per_loop, 4, NULL)) ;
        // GxB_Matrix_fprint(M, "M", GxB_SHORT, stdout);

        // // Remove rows with less than 3 entries
        // GRB_TRY (GrB_mxm (
        //     M_xor, NULL, NULL, bxor_first, M, hash_vals, NULL
        // )) ; 
        // GRB_TRY (GrB_Matrix_reduce_Monoid(M_outdeg, NULL, NULL, 
        //     GrB_MIN_MONOID_UINT64, M_xor, NULL)) ;
        // GRB_TRY (GrB_Vector_select_UINT8(
        //     M_outdeg, NULL, NULL, GrB_VALUENE_UINT8, M_outdeg, (uint8_t) 0, NULL
        // )) ;
        // GRB_TRY (GrB_Vector_nvals(&n_keep, M_outdeg));
        GrB_Index arr_size, junk_size;
        // if(n_keep == swaps_per_loop)
        // {
        //     M_fours = M;
        //     pairs_4s = P;
        // }
        // else // Remove bad swaps
        // {
        //     GRB_TRY (GxB_Vector_unpack_CSC(
        //         M_outdeg, &arr_keep, &junk, &arr_size, &junk_size, &iso,
        //         &n_keep, NULL, NULL)) ;
        //     LG_TRY (LAGraph_Free(&junk, msg));
            
        //     GRB_TRY (GrB_Matrix_new(&M_fours, GrB_UINT8, n_keep, 2)) ;
        //     GRB_TRY (GrB_Matrix_new(&pairs_4s, GrB_UINT8, n_keep, 2)) ; 

        //     GRB_TRY (GrB_Matrix_extract(
        //         M_fours, NULL, NULL, M, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
        //     GRB_TRY (GrB_Matrix_extract(
        //         pairs_4s, NULL, NULL, P, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
        //     GRB_TRY (GrB_free(&M));
        //     GRB_TRY (GrB_free(&P)); 
        //     LG_TRY (LAGraph_Free((void **)&arr_keep, msg));
        // }


        GRB_TRY (GrB_Matrix_new(&new_hashed_edges, GrB_UINT64, swaps_per_loop, 2)) ;
        GRB_TRY (GrB_Vector_new(&hashed_edges, GrB_UINT64, e)) ;

        //This Scalar is used in the hash funtion and should be something near 
        // 2^60 and odd.
        GRB_TRY (GrB_Matrix_assign_UINT64(
            hash_s, hash_s, NULL, 0xFB21C651E98DF25ULL, 
            GrB_ALL, 0, GrB_ALL, 0, GrB_DESC_S
        )) ;
        GRB_TRY (GrB_Vector_assign_UINT64(
            dense_hash, NULL, NULL, 0xFB21C651E98DF25ULL, GrB_ALL, 0, NULL
        )) ;
        GxB_Matrix_fprint(hash_s, "hash_s", GxB_SHORT, stdout);
        GRB_TRY(GrB_mxm(
            new_hashed_edges, NULL, NULL, hash_edges, 
            M, hash_s, NULL
        ));
        GRB_TRY(GrB_mxv(
            hashed_edges, NULL, NULL, hash_edges, 
            E, dense_hash, NULL
        ));
        // GxB_Vector_fprint(hashed_edges,"Hashes",GxB_SHORT, stdout);
        // GxB_Matrix_fprint(M_fours,"M",GxB_COMPLETE, stdout);
        // GxB_Matrix_fprint(new_hashed_edges,"Hashes",GxB_COMPLETE, stdout);
        

        // I will unpack and then reconstruct with hash as index.
        GRB_TRY(GxB_Matrix_unpack_FullR(
            new_hashed_edges, (void **) &hash_vals_new, &junk_size, &iso, NULL
        )) ;
        GRB_TRY(GxB_Vector_unpack_Full(
            hashed_edges, (void **) &hash_vals, &junk_size, &iso, NULL
        )) ;

        GRB_TRY (GrB_Vector_new(&exists, GrB_UINT64, 1ULL << 60)) ;
        GRB_TRY (GrB_Matrix_new(
            &buckets, GrB_UINT64, 1ULL << 60, swaps_per_loop)) ;
        
        // Build hash buckets
        GRB_TRY(GxB_Vector_build_Scalar(
            exists, hash_vals, one64, e
        )) ;
        // GxB_Vector_fprint(exists,"exists",GxB_SHORT, stdout);
        GRB_TRY(GxB_Matrix_build_Scalar(
            buckets, hash_vals_new, half_ramp, one64, swaps_per_loop * 2
        )) ;

        GRB_TRY(GrB_Matrix_reduce_Monoid(
            exists, NULL, GrB_PLUS_UINT64, GrB_PLUS_MONOID_UINT64, buckets, NULL
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
        
        GxB_Vector_fprint(r_exists,"r_exists",GxB_SHORT, stdout);
        GRB_TRY (GrB_Vector_assign_INT8(
            r_exists, r_exists, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC)) ;
        GRB_TRY (GxB_Vector_unpack_CSC(
                 r_exists, &arr_keep, &junk, &arr_size, &junk_size, &iso, &n_keep,
                 NULL,NULL));
        GRB_TRY(GrB_Matrix_new(&M_fours, GrB_UINT64, n_keep, 4)) ;
        //increase width so it can be used as a mask.
        GRB_TRY (GrB_Matrix_extract(
            M_fours, NULL, NULL, M, arr_keep, n_keep, 
            GrB_ALL, 0, NULL)) ;
        LG_TRY (LAGraph_Free(&junk, msg)) ;
        GRB_TRY (GrB_mxm(M_fours, NULL, NULL, GxB_ANY_FIRST_UINT64, M_fours, swap_p, NULL)) ;
        GRB_TRY (GrB_assign(M, NULL, NULL, M_fours, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GxB_Matrix_reshape(M, false, swaps_per_loop * 2, 2, NULL));
        
        GrB_Index I_range [3] ;
        I_range [GxB_BEGIN] = 0 ;
        I_range [GxB_END ] = swaps_per_loop * 2 - 1 ;
        GRB_TRY (GrB_assign(
            E, NULL, NULL, M, I_range, GxB_RANGE, GrB_ALL, 0, NULL
        )) ;
        GxB_Matrix_fprint(M, "M", GxB_SHORT, stdout);
        GxB_Matrix_fprint(E, "E", GxB_SHORT, stdout);
        
        // if(badcount == 0)
        // {
        //     E_arranged[1] = M_fours;
        //     pairs_new = pairs_4s;
        //     GRB_TRY (GrB_Matrix_new(E_arranged + 2, GrB_UINT8, n_keep, 2)) ;
        // }
        // else
        // {
        //     GRB_TRY(GrB_Vector_new(&r_exists, GrB_UINT8, n_keep)) ;
        //     // Search through array for bad swaps.
        //     GRB_TRY(GrB_vxm(
        //         r_exists, NULL, NULL, LAGraph_any_one_uint8, exists, buckets, NULL
        //     )) ;
            
        //     // Get compliment vector
        //     GRB_TRY (GrB_Vector_assign_INT64(
        //         r_exists, r_exists, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC)) ;
        //     GRB_TRY (GxB_Vector_unpack_CSC(
        //         r_exists, &arr_keep, &junk, &arr_size, &junk_size, &iso, &n_keep,
        //         NULL, NULL)) ;
        //     LG_TRY (LAGraph_Free(&junk, msg)) ;
            
        //     GRB_TRY (GrB_Matrix_new(E_arranged + 1, GrB_UINT8, n_keep, 2)) ;
        //     GRB_TRY (GrB_Matrix_new(E_arranged + 2, GrB_UINT8, n_keep, 2)) ;
        //     GRB_TRY (GrB_Matrix_new(&pairs_new, GrB_UINT8, n_keep, e)) ;

        //     GRB_TRY (GrB_Matrix_extract(
        //         E_arranged[1], NULL, NULL, M_fours, arr_keep, n_keep, 
        //         GrB_ALL, 0, NULL)) ;
        //     GRB_TRY (GrB_Matrix_extract(
        //         pairs_new, NULL, NULL, pairs_4s, arr_keep, n_keep, 
        //         GrB_ALL, 0, NULL)) ;
        //     LG_TRY (LAGraph_Free((void **)&arr_keep, msg));
        // }

        //Split the swapped edges into their correct places.
        // GRB_TRY (GrB_select(
        //     E_arranged[2], NULL, NULL, first_bit, 
        //     E_arranged[1], (uint8_t) 0, NULL
        // )) ;
        // GRB_TRY (GrB_select(
        //     E_arranged[1], NULL, NULL, first_bit, 
        //     E_arranged[1], (uint8_t) 1, GrB_DESC_R
        // )) ;

        // GRB_TRY (GrB_Matrix_apply_BinaryOp2nd_UINT8(
        //     E_arranged[1], NULL, NULL, GrB_DIV_UINT8, E_arranged[1], 2, NULL
        // )) ;
        // GRB_TRY (GrB_Matrix_apply_BinaryOp2nd_UINT8(
        //     E_arranged[2], NULL, NULL, GrB_DIV_UINT8, E_arranged[2], 2, NULL
        // )) ;
        // GxB_Matrix_fprint(E_arranged[1],"news",GxB_COMPLETE, stdout);
        // GxB_Matrix_fprint(E_arranged[2],"news",GxB_COMPLETE, stdout);
        // GrB_Index n_old;
        // GRB_TRY (GrB_Matrix_reduce_Monoid(
        //     r_pairs, NULL, NULL, GxB_ANY_UINT8_MONOID, pairs_new, GrB_DESC_RT0));

        // GRB_TRY (GrB_Vector_assign_UINT8(
        //     r_pairs, r_pairs, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC)) ;
        
        // GRB_TRY (GxB_Vector_unpack_CSC(
        //     r_pairs, &arr_keep, &junk, &arr_size, &junk_size, &iso, &n_old,
        //     NULL, NULL)) ;
        // LG_TRY (LAGraph_Free(&junk, msg)) ;
        
        // LG_ASSERT(n_old ==  e - 2 * n_keep, GrB_NOT_IMPLEMENTED) ;

        // // old edges that are being kept.
        // GRB_TRY (GrB_Matrix_new(E_arranged, GrB_UINT8, n_old, n)) ;
        // GRB_TRY (GrB_Matrix_extract(E_arranged[0], NULL, NULL, E, 
        //     arr_keep, n_old, GrB_ALL, n, NULL)) ;
        // LG_TRY (LAGraph_Free((void **)&arr_keep, msg));
        
        // E = Concat(E_prime, M_1, M_2)
        // where E_prime contains no edges in the indegree of pairs after it is 
        // pruned
        // GRB_TRY (GxB_Matrix_concat(E, E_arranged, 3, 1, NULL)) ;
        
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
        E, (void **)&row_indices, &ind_size, &iso, NULL)) ;
    GRB_TRY (GxB_Matrix_build_Scalar(*A_new, row_indices, row_indices + e, one8, e));
    GRB_TRY (GrB_eWiseAdd(*A_new, NULL, NULL, GrB_PLUS_UINT8, *A_new,*A_new, GrB_DESC_T0)) ;
    GxB_Matrix_fprint(*A_new,"A_new",GxB_SHORT, stdout);
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
