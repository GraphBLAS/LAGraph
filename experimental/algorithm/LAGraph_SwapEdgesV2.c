//------------------------------------------------------------------------------
// LAGraph_SwapEdges: randomly swaps edges in a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel Gomez, Texas A&M University

//------------------------------------------------------------------------------

// References:

// R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, and U. Alon, “On the 
// uniform generation of random graphs with prescribed degree sequences,” 2004.

#define FREE_LOOP                               \
{                                               \
    GrB_free (&M) ;                             \
    GrB_free(&dup_swaps_v);                     \
    GrB_free(&new_hashed_edges);                \
    GrB_free(&hashed_edges);                    \
    GrB_free(&edge_perm) ;                      \
}

#define LG_FREE_WORK                            \
{                                               \
    /* free any workspace used here */          \
    GrB_free (&E) ;                             \
    GrB_free (&A_tril) ;                        \
    GrB_free (&random_v) ;                      \
    GrB_free (&r_permute) ;                     \
    GrB_free (&ramp_v) ;                        \
    GrB_free (&E_temp) ;                        \
    GrB_free (&r_60) ;                          \
    GrB_free(&exists);                          \
    GrB_free (&E_vec) ;                         \
    GrB_free (&swap_pair) ;                     \
    GrB_free (&hash_seed_e) ;                   \
    GrB_free (&add_term_biop) ;                 \
    GrB_free (&add_term_monoid) ;               \
    GrB_free (&plus_term_one) ;                 \
    GrB_free (&second_edge) ;                   \
    GrB_free (&second_bool_edge) ;              \
    GrB_free (&second_edge_monoid) ;            \
    GrB_free (&second_second_edge) ;            \
    GrB_free (&lg_shiftland) ;                  \
    GrB_free (&lg_edge) ;                       \
    GrB_free (&lg_swap) ;                       \
    GrB_free (&con) ;                           \
    GrB_free (&E_con) ;                         \
    GrB_free (&one8) ;                          \
    GrB_free (&x) ;                             \
    LAGraph_Free((void**)&indices, msg) ;       \
    LAGraph_Free((void **) &dup_swaps, NULL);   \
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
#if GxB_IMPLEMENTATION >= GxB_VERSION (10,0,0)
void shift_and 
    (uint16_t *z, const uint16_t *x)
    {
        (*z) = (*x) & ((*x) << 8);
        (*z) |= (*z) >> 8;
    }
#define SHIFT_AND                                                               \
"void shift_and                                                              \n"\
"   (uint16_t *z, const uint16_t *x)                                         \n"\
"   {                                                                        \n"\
"       (*z) = (*x) & ((*x) << 8);                                           \n"\
"       (*z) |= (*z) >> 8;                                                   \n"\
"   }"

typedef struct {
    uint64_t a; 
    uint64_t b;
} edge_type;
#define EDGE_TYPE                                                               \
"typedef struct { uint64_t a; uint64_t b; } edge_type;"

typedef struct {
    uint64_t a; 
    uint64_t b;
    uint64_t c; 
    uint64_t d;
} swap_type;
#define SWAP_TYPE                                                               \
"typedef struct {                                                            \n"\
"   uint64_t a; uint64_t b; uint64_t c; uint64_t d;                         \n" \
"}swap_type;"

void swap_bc
(swap_type *z, const swap_type *x, GrB_Index I, GrB_Index J, const bool *y)
{
    memcpy(z, x, sizeof(*z)) ; //unnessesary when aliassed but done for safety.
    if(z->a == z->c || z->b == z->c || z->a == z->d || z->b == z->d ) return;
    if(I & 1)
    {
        uint64_t temp = z->d;
        z->d = z->b;
        z->b = temp; 
    }
    else
    {
        uint64_t temp = z->c;
        z->c = z->b;
        z->b = temp; 
    }    
}
#define SWAP_BC                                                                 \
"void swap_bc                                                                   \n"\
"(swap_type *z, const swap_type *x, GrB_Index I, GrB_Index J, const bool *y)    \n"\
"{                                                                              \n"\
"    memcpy(z, x, sizeof(*z)) ; //unnessesary when aliassed but done for safety. \n"\
"    if(z->a == z->c || z->b == z->c || z->a == z->d || z->b == z->d ) return;  \n"\
"   if(I & 1)                                                                   \n"\
"    {                                                                          \n"\
"        uint64_t temp = z->d;                                                  \n"\
"        z->d = z->b;                                                           \n"\
"        z->b = temp;                                                           \n"\
"    }                                                                          \n"\
"    else                                                                       \n"\
"    {                                                                          \n"\
"        uint64_t temp = z->c;                                                  \n"\
"        z->c = z->b;                                                           \n"\
"        z->b = temp;                                                           \n"\
"    }                                                                          \n"\
"}"


//Simply making a cantor pairing then masking.
void hash_edge 
(uint64_t *z, const edge_type *x, const uint64_t *mask)
{
    (*z) = (((x->a + x->b + 1) * (x->a + x->b)) / 2) & (*mask) ;
    (*z) += (x->a < x->b)? x->a: x->b;
    (*z) &= (*mask);
}
#define HASH_EDGE                                                                \
"void hash_edge                                                               \n"\
"(uint64_t *z, const edge_type *x, const uint64_t *mask)                      \n"\
"{                                                                            \n"\
"    (*z) = (((x->a + x->b + 1) * (x->a + x->b)) / 2) & (*mask) ;             \n"\
"    (*z) += (x->a < x->b)? x->a: x->b;                                       \n"\
"    (*z) &= (*mask);                                                         \n"\
"}"

void add_term
    (uint8_t *z, const uint8_t *x, const uint8_t *y)
{
    (*z) = (*x) | (*y) + ((uint8_t)1 & (*x) & (*y)) ;
}
#define ADD_TERM                                                               \
"void add_term                                                                \n"\
"(uint8_t *z, const uint8_t *x, const uint8_t *y)                             \n"\
"{                                                                            \n"\
"    (*z) = (*x) | (*y) + ((uint8_t)1 & (*x) & (*y)) ;                         \n"\
"}"
void edge2
    (edge_type *z, const void *x, const edge_type *y)
{
    //if(y->a == 0 && y->b == 0) return;
    z->a = y->a;
    z->b = y->b;
}
#define EDGE2                                                                   \
"void edge2                                                                   \n"\
"(edge_type *z, const void *x, const edge_type *y)                            \n"\
"{                                                                            \n"\
"    //if(y->a == 0 && y->b == 0) return;                                     \n"\
"    z->a = y->a;                                                             \n"\
"    z->b = y->b;                                                             \n"\
"}"
// void edge2b
//     (edge_type *z, const void *x, const edge_type *y)
// {
//     z->a = y->a;
//     z->b = y->b;
// }
// #define EDGE2B                                                                   \
// "void edge2b                                                                   \n"\
// "(edge_type *z, const void *x, const edge_type *y)                            \n"\
// "{                                                                            \n"\
// "    z->a = y->a;                                                             \n"\
// "    z->b = y->b;                                                             \n"\
// "}"

int LAGraph_SwapEdgesV2
(
    // output
    GrB_Matrix *A_new, //The adjacency matrix of G with edges randomly swapped
    // input: not modified
    LAGraph_Graph G,
    GrB_Index Q, // Swaps per edge
    char *msg
)
{
    // TODO: Remove all instances of GxB_(un)pack!
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL; // n x n Adjacency Matrix 

    // e x 2 with entries corresponding to verticies of an edge
    GrB_Matrix E = NULL;
    // e x 1 vector, each entry is an edge.
    GrB_Vector E_vec = NULL, E_temp = NULL; 
    GxB_Container E_con = NULL;

    // swaps x 4
    // Each row contains 4 entries corresponding to the verticies 
    // that are involved in the swap.
    GrB_Vector M = NULL;
    GxB_Container con = NULL; 

    // n = |V| e = |E|
    GrB_Index n = 0, e = 0;

    // n x n 
    // Lower triangle of adjacency matrix
    GrB_Matrix A_tril = NULL ;

    // e x 1 random vectors
    GrB_Vector random_v = NULL, r_permute = NULL;

    // indicies for A
    GrB_Index *indices = NULL;

    GrB_Vector ramp_v = NULL;

    // edge permutation
    GrB_Vector edge_perm = NULL;
    bool iso = false;

    // Number of values kept in each phase
    GrB_Index n_keep = 0, M_nvals = 0, dup_arr_size = 0;

    // swaps x 2 matrix which holds the hashes of each planned edge. 
    GrB_Vector new_hashed_edges = NULL;

    // e holds hashes of old edges
    GrB_Vector hashed_edges = NULL;

    // 2^60 holds the buckets in which hashes collided.
    GrB_Vector exists = NULL; 
    
    GrB_UnaryOp lg_shiftland = NULL;

    //  b1 <---> a2 or b1 <---> b2
    GrB_IndexUnaryOp swap_pair = NULL;
    
    // z = h_y(x)
    GrB_BinaryOp hash_seed_e = NULL;

    // This monoid has only been designed for inputs in {0,1,2}, other behavior 
    // is undefined.
    // (0,x) -> x, (1,1) -> 2, (2,x) -> 2 (and commutative)
    // Aka z = min(2, x + y)
    GrB_BinaryOp add_term_biop = NULL;
    GrB_Monoid add_term_monoid = NULL;
    GrB_Semiring plus_term_one = NULL;

    // z = y
    GrB_BinaryOp second_edge = NULL;
    GrB_BinaryOp second_bool_edge = NULL;
    GrB_Monoid second_edge_monoid = NULL;
    GrB_Semiring second_second_edge = NULL;

    // Toople types
    GrB_Type lg_edge = NULL, lg_swap = NULL;
    // Unload types
    GrB_Type M_type = NULL, E_type = NULL;
    int M_hand = 0, E_hand = 0;

    int16_t *dup_swaps = NULL;
    GrB_Vector dup_swaps_v = NULL;
    // BOOL swaps * 2 vector that holds false if an edge in the swap did not "work"

    GrB_Vector sort_h = NULL;
    GrB_Vector r_60 = NULL;

    // Constants ---------------------------------------------------------------
    GrB_Vector x = NULL;

    GrB_Scalar one8 = NULL;

    GrB_Index ind_size = 0;
    
    //--------------------------------------------------------------------------
    // Check inputs TODO
    //--------------------------------------------------------------------------
    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED,
        LAGRAPH_INVALID_GRAPH, 
        "G must be undirected"
    ) ;
    // char type[LAGRAPH_MAX_NAME_LEN];
    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, 
        "G->nself_edges must be zero") ;
    LG_ASSERT (A_new != NULL, GrB_NULL_POINTER) ;
    // GRB_TRY (GxB_Matrix_type(&E_type, A)) ;
    // LG_ASSERT_MSG (E_type == GrB_BOOL || , LAGRAPH_INVALID_GRAPH, 
    //     "A must be structural") ;

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;  

    // Types
    GRB_TRY (GxB_Type_new(
        &lg_edge, sizeof(edge_type), "edge_type", EDGE_TYPE)) ;
    GRB_TRY (GxB_Type_new(
        &lg_swap, sizeof(swap_type), "swap_type", SWAP_TYPE)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new (&A_tril, GrB_BOOL, n, n)) ;


    // Extract lower triangular edges.
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_nvals(&e, A_tril)) ;

    
    GRB_TRY (GrB_Matrix_new(&E, GrB_UINT64, e, 2)) ;
    GRB_TRY (GrB_Vector_new(&E_vec, lg_edge, e)) ;
        
    //Init Operators -----------------------------------------------------------
    GRB_TRY (GxB_UnaryOp_new (
        &lg_shiftland, (GxB_unary_function) (&shift_and),
        GrB_UINT16, GrB_UINT16, "shift_and", SHIFT_AND
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &hash_seed_e, (GxB_binary_function) (&hash_edge),
        GrB_UINT64, lg_edge, GrB_UINT64, "hash_edge", HASH_EDGE
    )) ;
    GRB_TRY (GxB_IndexUnaryOp_new (
        &swap_pair, (GxB_index_unary_function) (&swap_bc),
        lg_swap, lg_swap, GrB_BOOL, "swap_bc", SWAP_BC
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &add_term_biop, (GxB_binary_function) (&add_term), 
        GrB_UINT8, GrB_UINT8, GrB_UINT8, "add_term", ADD_TERM
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &second_edge, (GxB_binary_function) (&edge2), 
        lg_edge, lg_edge, lg_edge, "edge2", EDGE2
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &second_bool_edge, (GxB_binary_function) (&edge2), 
        lg_edge, GrB_BOOL, lg_edge, "edge2", EDGE2
    )) ;

    GRB_TRY (GxB_Monoid_terminal_new_UINT8(
        &add_term_monoid, add_term_biop, (uint8_t) 0, (uint8_t) 2
    )) ;

    // This isn't actually a monoid but since it's never applied as one it 
    // shouldn't be a problem? FORESHADOWING
    edge_type iden_second = {0,0};
    GRB_TRY (GrB_Monoid_new_UDT(
        &second_edge_monoid, second_edge, (void *) &iden_second
    )) ;

    // Now working with the built-in ONEB binary op
    GRB_TRY(GrB_Semiring_new(
        &plus_term_one, add_term_monoid, GrB_ONEB_UINT8
    )) ;
    GRB_TRY(GrB_Semiring_new(
        &second_second_edge, second_edge_monoid, second_bool_edge
    )) ;
    // count swaps 
    GrB_Index num_swaps = 0, num_attempts = 0, swaps_per_loop = e / 3 ;

    // Make E Matrix -----------------------------------------------------------
    LG_TRY (LAGraph_Malloc (
        (void**)(&indices), 2ull * e, sizeof(GrB_Index), msg)) ;
    GRB_TRY (
        GrB_Matrix_extractTuples_BOOL (indices, indices + e, NULL, &e, A_tril)
        ) ;
    GRB_TRY (GxB_Matrix_pack_FullC (
        E, (void **)&indices, 2ull * e * sizeof(GrB_Index), false, NULL
    )) ;
    GRB_TRY (GrB_set(E, GrB_ROWMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
    int shift_e = __builtin_clzl(e);
    uint64_t ehash_size = (1ull << (67-shift_e)) ;
    // if(ehash_size > 6*e) ehash_size/=2;
    printf("Hash Size: %ld", ehash_size);
    GRB_TRY (GrB_Vector_new(&exists, GrB_UINT8, ehash_size)) ;
    
    // GxB_Matrix_fprint(E, "E", GxB_SHORT, stdout);
    // Init Ramps --------------------------------------------------------------
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, e + 1)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (ramp_v, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    GrB_Index ramp_size;

    // Init Constants ----------------------------------------------------------
    GRB_TRY (GrB_Scalar_new (&one8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (one8, 1)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_BOOL, e)) ;

    // Make Random -------------------------------------------------------------
    GRB_TRY (GrB_Vector_new(&random_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_60, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, 1ull << (64-shift_e))) ;
    GRB_TRY(GrB_set (r_permute, GxB_BITMAP, GxB_SPARSITY_CONTROL)) ;
    GRB_TRY(GrB_set (exists, GxB_BITMAP | GxB_FULL, GxB_SPARSITY_CONTROL)) ;
    // GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        random_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;
    //TODO: Change seed
    LG_TRY(
        LAGraph_Random_Seed(random_v, 1548945616ul, msg)) ;
    GRB_TRY (GxB_Container_new(&E_con)) ;
    GRB_TRY (GxB_unload_Matrix_into_Container(E, E_con, NULL));
    GRB_TRY (GxB_Vector_unload(
        E_con->x, (void **) &indices, &E_type, &e, &ind_size, &E_hand, NULL));
    e /= 2;
    GRB_TRY (GxB_Vector_load(
        E_vec, (void **) &indices, lg_edge, e, ind_size, E_hand, NULL));
    printf("Entering loop, Good Luck:\n") ;
    while(num_swaps < e * Q && num_attempts < e * Q * 5)
    {
        GrB_Index perm_size, arr_size, junk_size;
        // Coming into the loop: 
        // E must be the incidence matrix of the new graph. W/o self edges nor 
        // parallel edges. Each row must have exactly two distinct values.
        // random_v has a radom dense vector.
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
            r_60, NULL, NULL, GxB_BSHIFT_UINT64, random_v, -(shift_e), NULL
        )) ;
        GRB_TRY (GrB_Vector_clear(x)) ;
        GRB_TRY (GrB_Vector_resize(x, e)) ;
        GRB_TRY (GrB_Vector_assign_BOOL(
            x, NULL, NULL, true, GrB_ALL, 0, NULL)) ;
        LG_TRY (LAGraph_FastAssign(
            r_permute, NULL, NULL, r_60, x, ramp_v, GxB_ANY_FIRSTJ_INT64,
            NULL, msg)) ;
        
        GrB_Index edges_permed = 0;
        GRB_TRY (GrB_Vector_nvals(&edges_permed, r_permute)) ;
        GRB_TRY (GrB_Vector_new(&edge_perm, GrB_BOOL, edges_permed)) ;
        GRB_TRY (GrB_Vector_extractTuples(NULL, edge_perm, r_permute, NULL)) ; 
        swaps_per_loop = LAGRAPH_MIN(swaps_per_loop, edges_permed / 2) ;

        // Chose only the edges we need from  edge_perm
        GRB_TRY (GrB_Vector_resize(edge_perm, swaps_per_loop * 2)) ;
        
        // Get the desired edges from the E_vec array.
        GRB_TRY (GrB_Vector_new(&M, lg_edge, swaps_per_loop * 2)) ;
        GRB_TRY (GxB_Vector_extract_Vector(
            M, NULL, NULL, E_vec, edge_perm, NULL
        )) ; 

        // Make the swaps via the swap_pair unary op.
        GRB_TRY (GxB_Vector_unload(
            M, (void **) &indices, &M_type, &M_nvals, &ind_size, &M_hand,
            NULL
        )) ;
        GRB_TRY (GxB_Vector_load(
            M, (void **) &indices, lg_swap, M_nvals / 2, ind_size, M_hand, NULL
        )) ;
        GRB_TRY (GrB_Vector_apply_IndexOp_BOOL(
            M, NULL, NULL, swap_pair, M, false, NULL)) ;
        GRB_TRY (GxB_Vector_unload(
            M, (void **) &indices, &M_type, &M_nvals, &ind_size, &M_hand,
            NULL
        )) ;
        GRB_TRY (GxB_Vector_load(
            M, (void **) &indices, lg_edge, M_nvals * 2, ind_size, M_hand, 
            NULL
        )) ;

        // Hash Edges ----------------------------------------------------------
        GRB_TRY (GrB_Vector_new(
            &new_hashed_edges, GrB_UINT64, swaps_per_loop * 2)) ;
        GRB_TRY (GrB_Vector_new(&hashed_edges, GrB_UINT64, e)) ;

        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
            new_hashed_edges, NULL, NULL, hash_seed_e, M, 
            ehash_size - 1ll, NULL
        )) ;//0xFB21C651E98DF25ULL
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
            hashed_edges, NULL, NULL, hash_seed_e, E_vec, 
            ehash_size - 1ll, NULL
        )) ;        

        //----------------------------------------------------------------------
        // Build Hash Buckets
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Vector_new(&dup_swaps_v, GrB_INT8, swaps_per_loop * 2)) ;
        GRB_TRY (GrB_set(dup_swaps_v, GxB_BITMAP, GxB_SPARSITY_CONTROL)) ;
        
        GRB_TRY (GrB_Vector_clear(x)) ;
        GRB_TRY (GrB_Vector_resize(x, e)) ;
        GRB_TRY (GrB_Vector_assign_BOOL(
            x, NULL, NULL, true, GrB_ALL, 0, NULL)) ;
        // place a one in any bucket that coresponds to an edge currently in E
        LG_TRY (LAGraph_FastAssign(
            exists, NULL, NULL, hashed_edges, x, ramp_v, GxB_ANY_PAIR_UINT8,
            NULL, msg
        )) ;
        
        GRB_TRY (GrB_Vector_clear(x)) ;
        GRB_TRY (GrB_Vector_resize(x, swaps_per_loop * 2)) ;
        GRB_TRY (GrB_Vector_assign_BOOL(
            x, NULL, NULL, true, GrB_ALL, 0, NULL)) ;

        // exists cannot possibly be full at this point.
        // Fill out exists in O(1) time.
        GRB_TRY (GxB_Container_new(&con)) ;
        GRB_TRY (GxB_unload_Vector_into_Container(exists, con, NULL)) ;
        // Sanity check
        LG_ASSERT (con->format == GxB_BITMAP, GrB_INVALID_VALUE) ;
        GrB_free(&exists);
        exists = con->b;
        con->b = NULL;
        GRB_TRY (GrB_free(&con)) ;
        GRB_TRY(GrB_set (exists, GxB_BITMAP | GxB_FULL, GxB_SPARSITY_CONTROL)) ;

        // "Count" all of the edges that fit into each bucket. Stop counting at 
        // 2 since we will have to throw that whole bucket away anyway.
        LG_TRY (LAGraph_FastAssign(
            exists, NULL, add_term_biop, new_hashed_edges, x, ramp_v, 
            plus_term_one, NULL, msg
        )) ;

        // Select buckets with only one corresponding value
        GRB_TRY (GrB_Vector_select_UINT8(
            exists, NULL, NULL, GrB_VALUEEQ_UINT8, exists, (uint8_t) 1,
            NULL
        )) ;

        // Find each hashed edge's bucket, dup_swaps_v is 1 if exists[edge] = 1
        LG_TRY (LAGraph_FastAssign(
            dup_swaps_v, NULL, NULL, new_hashed_edges, exists, ramp_v, 
            GxB_ANY_PAIR_INT8, GrB_DESC_T0, msg
        )) ;

        // Fill out dup_swaps_v in O(1) time.
        GRB_TRY (GxB_Container_new(&con)) ;
        GRB_TRY (GxB_unload_Vector_into_Container(dup_swaps_v, con, NULL)) ;
        n_keep = con->nvals;
        GRB_TRY (GrB_free(&dup_swaps_v));
        dup_swaps_v = con->b;
        con->b = NULL;
        GRB_TRY (GrB_free(&con)) ;
        GRB_TRY (GxB_Vector_unload(
            dup_swaps_v, (void **) &dup_swaps, &M_type, &M_nvals, &dup_arr_size, 
            &M_hand, NULL)) ;
        GRB_TRY (GxB_Vector_load(
            dup_swaps_v, (void **) &dup_swaps, GrB_INT16, M_nvals / 2, 
            dup_arr_size, M_hand, NULL)) ;
        GRB_TRY (GrB_apply(
            dup_swaps_v, NULL, NULL, lg_shiftland, dup_swaps_v, NULL)) ;
        GRB_TRY (GxB_Vector_unload(
            dup_swaps_v, (void **) &dup_swaps, &M_type, &M_nvals, &dup_arr_size, 
            &M_hand, NULL)) ;
        GRB_TRY (GxB_Vector_load(
            dup_swaps_v, (void **) &dup_swaps, GrB_INT8, M_nvals * 2, 
            dup_arr_size, M_hand, NULL)) ;
        
        GRB_TRY (GrB_Vector_clear(exists)) ;
        // ---------------------------------------------------------------------
        // Place Good Swaps back into E_vec
        // ---------------------------------------------------------------------

        GRB_TRY (GxB_Container_new(&con)) ;
        GRB_TRY (GxB_unload_Vector_into_Container(M, con, NULL)) ;
        GRB_TRY (GrB_free(&(con->b))) ;
        con->b = dup_swaps_v;
        con->nvals = n_keep;
        con->format = GxB_BITMAP;
        dup_swaps_v = NULL;
        GRB_TRY (GxB_load_Vector_from_Container(M, con, NULL)) ;
        GRB_TRY (GrB_free(&con)) ;
        GRB_TRY (GrB_Vector_new(&E_temp, lg_edge, e)) ;

        // Want to place edges of M into E_vec inplace BUT:
        // EVIL segfault
        // GRB_TRY (LAGraph_FastAssign(
        //     E_vec, NULL, second_edge, edge_perm, M, ramp_v, 
        //     second_second_edge, NULL, msg)) ;

        // Temporary fix for the above:
        GRB_TRY (LAGraph_FastAssign(
            E_temp, NULL, NULL, edge_perm, M, ramp_v, 
            second_second_edge, NULL, msg)) ;
        GRB_TRY (GrB_eWiseAdd(
            E_vec, NULL, NULL, second_edge, E_vec, E_temp, NULL)) ;
        GrB_free (&E_temp) ;

        n_keep /= 2;

        FREE_LOOP ; // Free Matricies that have to be rebuilt

        // Adjust number of swaps to do next.
        num_attempts += swaps_per_loop ;
        num_swaps += n_keep ;
        swaps_per_loop = n_keep * 3 ;
        swaps_per_loop = LAGRAPH_MAX(swaps_per_loop, 16) ;
        swaps_per_loop = LAGRAPH_MIN(swaps_per_loop, e / 3) ;

        LG_TRY (LAGraph_Random_Next(random_v, msg)) ;
        printf("#####Made %ld swaps. Total %ld out of %ld. Attempting %ld swaps next.#####\n\n", n_keep, num_swaps, e * Q, swaps_per_loop) ;
    } 
    GRB_TRY (GxB_Vector_unload(
        E_vec, (void **) &indices, &lg_edge, &e, &ind_size, &E_hand, NULL));
    GRB_TRY (GxB_Vector_load(
        E_con->x, (void **) &indices, E_type, e * 2, ind_size, E_hand, NULL));
    GRB_TRY (GxB_load_Matrix_from_Container(E,E_con, NULL)) ;
    GRB_TRY (GrB_set(E, GrB_COLMAJOR, GrB_STORAGE_ORIENTATION_HINT)) ;
    GRB_TRY (GxB_unload_Matrix_into_Container(E,E_con, NULL)) ;
    GRB_TRY (GxB_Vector_unload(
        E_con->x, (void **) &indices, &E_type, &e, &ind_size, &E_hand, NULL));
    e /= 2;

    // Build Output Matrix
    GRB_TRY(GrB_Matrix_new(A_new, GrB_BOOL, n, n)) ;
    GRB_TRY (GxB_Matrix_build_Scalar(*A_new, indices, indices + e, one8, e)) ;
    GRB_TRY (GrB_eWiseAdd(
        *A_new, NULL, NULL, GrB_LOR_MONOID_BOOL, *A_new,*A_new, GrB_DESC_T0
    )) ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
#endif