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
    GrB_free(&edge_perm) ;                      \
}

#define LG_FREE_WORK                            \
{                                               \
    /* free any workspace used here */          \
    GrB_free (&A_tril) ;                        \
    GrB_free (&Ai) ;                            \
    GrB_free (&Aj) ;                            \
    GrB_free (&random_v) ;                      \
    GrB_free (&r_permute) ;                     \
    GrB_free (&ramp_v) ;                        \
    GrB_free (&E_temp) ;                        \
    GrB_free (&r_60) ;                          \
    GrB_free (&exists) ;                        \
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
    GrB_free (&hashed_edges) ;                  \
    GrB_free (&con) ;                           \
    GrB_free (&one8) ;                          \
    GrB_free (&x) ;                             \
    GrB_free (&temp_i) ;                        \
    LAGraph_Free ((void ** )&indices, NULL) ;   \
    LAGraph_Free ((void **) &dup_swaps, NULL);  \
    FREE_LOOP ;                                 \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (&A_new) ;                     \
    LAGraph_Delete(G_new, NULL) ;           \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"
void LG_SE_shift_and 
    (uint16_t *z, const uint16_t *x)
    {
        (*z) = (*x) & ((*x) << 8);
        (*z) |= (*z) >> 8;
    }
#define SHIFT_AND                                                               \
"void LG_SE_shift_and                                                        \n"\
"   (uint16_t *z, const uint16_t *x)                                         \n"\
"   {                                                                        \n"\
"       (*z) = (*x) & ((*x) << 8);                                           \n"\
"       (*z) |= (*z) >> 8;                                                   \n"\
"   }"

typedef struct {
    uint64_t a; 
    uint64_t b;
} LG_SE_edge_type64;
#define EDGE_TYPE64                                                             \
"typedef struct { uint64_t a; uint64_t b; } LG_SE_edge_type64;"

typedef struct {
    uint64_t a; 
    uint64_t b;
    uint64_t c; 
    uint64_t d;
} LG_SE_swap_type64;
#define SWAP_TYPE64                                                             \
"typedef struct {                                                            \n"\
"   uint64_t a; uint64_t b; uint64_t c; uint64_t d;                          \n"\
"} LG_SE_swap_type64;"

typedef struct {
    uint32_t a; 
    uint32_t b;
} LG_SE_edge_type32;
#define EDGE_TYPE32                                                             \
"typedef struct { uint32_t a; uint32_t b; } LG_SE_edge_type32;"

typedef struct {
    uint32_t a; 
    uint32_t b;
    uint32_t c; 
    uint32_t d;
} LG_SE_swap_type32;
#define SWAP_TYPE32                                                             \
"typedef struct {                                                            \n"\
"   uint32_t a; uint32_t b; uint32_t c; uint32_t d;                          \n"\
"}LG_SE_swap_type32;"

void LG_SE_swap_bc64
(LG_SE_swap_type64 *z, const LG_SE_swap_type64 *x, GrB_Index I, GrB_Index J, const bool *y)
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
void LG_SE_swap_bc32
(LG_SE_swap_type32 *z, const LG_SE_swap_type32 *x, GrB_Index I, GrB_Index J, const bool *y)
{
    memcpy(z, x, sizeof(*z)) ; //unnessesary when aliassed but done for safety.
    if(z->a == z->c || z->b == z->c || z->a == z->d || z->b == z->d ) return;
    if(I & 1)
    {
        uint32_t temp = z->d;
        z->d = z->b;
        z->b = temp; 
    }
    else
    {
        uint32_t temp = z->c;
        z->c = z->b;
        z->b = temp; 
    }    
}
#define SWAP_BC64                                                                  \
"void LG_SE_swap_bc64                                                           \n"\
"(LG_SE_swap_type64 *z, const LG_SE_swap_type64 *x, GrB_Index I, GrB_Index J, const bool *y)\n"\
"{                                                                              \n"\
"    memcpy(z, x, sizeof(*z)) ; //unnessesary when aliassed but done for safety.\n"\
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
#define SWAP_BC32                                                                  \
"void LG_SE_swap_bc32                                                           \n"\
"(LG_SE_swap_type32 *z, const LG_SE_swap_type32 *x, GrB_Index I, GrB_Index J, const bool *y)\n"\
"{                                                                              \n"\
"    memcpy(z, x, sizeof(*z)) ; //unnessesary when aliassed but done for safety.\n"\
"    if(z->a == z->c || z->b == z->c || z->a == z->d || z->b == z->d ) return;  \n"\
"   if(I & 1)                                                                   \n"\
"    {                                                                          \n"\
"        uint32_t temp = z->d;                                                  \n"\
"        z->d = z->b;                                                           \n"\
"        z->b = temp;                                                           \n"\
"    }                                                                          \n"\
"    else                                                                       \n"\
"    {                                                                          \n"\
"        uint32_t temp = z->c;                                                  \n"\
"        z->c = z->b;                                                           \n"\
"        z->b = temp;                                                           \n"\
"    }                                                                          \n"\
"}"

// using xorshift, from https://en.wikipedia.org/wiki/Xorshift
// with a state of uint64_t, or xorshift64star.
void LG_SE_hash_edge64
(uint64_t *z, const LG_SE_edge_type64 *x, const uint64_t *mask)
{
    (*z) = x->a ^ x->b;
	(*z) ^= (*z) << 13;
	(*z) ^= (*z) >> 7;
    (*z) ^= (x->a < x->b)? x->a: x->b;
	(*z) ^= (*z) << 17;
    (*z) &= (*mask);
}
void LG_SE_hash_edge32
(uint64_t *z, const LG_SE_edge_type32 *x, const uint64_t *mask)
{
    (*z) = x->a ^ x->b;
	(*z) ^= (*z) << 13;
	(*z) ^= (*z) >> 7;
    (*z) ^= (uint64_t)((x->a < x->b)? x->a: x->b);
	(*z) ^= (*z) << 17;
    (*z) &= (*mask);
}
#define HASH_EDGE64                                                              \
"void LG_SE_hash_edge64                                                       \n"\
"(uint64_t *z, const LG_SE_edge_type64 *x, const uint64_t *mask)              \n"\
"{                                                                            \n"\
"   (*z) = x->a ^ x->b;                                                       \n"\
"	(*z) ^= (*z) << 13;                                                   \n"\
"	(*z) ^= (*z) >> 7;                                                    \n"\
"   (*z) ^= (uint64_t)((x->a < x->b)? x->a: x->b);                            \n"\
"	(*z) ^= (*z) << 17;                                                   \n"\
"   (*z) &= (*mask);                                                          \n"\
"}"
#define HASH_EDGE32                                                              \
"void LG_SE_hash_edge32                                                       \n"\
"(uint64_t *z, const LG_SE_edge_type32 *x, const uint64_t *mask)              \n"\
"{                                                                            \n"\
"   (*z) = x->a ^ x->b;                                                       \n"\
"	(*z) ^= (*z) << 13;                                                   \n"\
"	(*z) ^= (*z) >> 7;                                                    \n"\
"   (*z) ^= (uint64_t)((x->a < x->b)? x->a: x->b);                            \n"\
"	(*z) ^= (*z) << 17;                                                   \n"\
"   (*z) &= (*mask);                                                          \n"\
"}"

void LG_SE_add_term
    (int8_t *z, const int8_t *x, const int8_t *y)
{
    (*z) = (*x) | (*y) + ((int8_t)1 & (*x) & (*y)) ;
}
#define ADD_TERM                                                                 \
"void LG_SE_add_term                                                          \n"\
"(int8_t *z, const int8_t *x, const int8_t *y)                                \n"\
"{                                                                            \n"\
"    (*z) = (*x) | (*y) + ((int8_t)1 & (*x) & (*y)) ;                         \n"\
"}"

void LG_SE_edge2nd64_bool
    (LG_SE_edge_type64 *z, const bool *x, const LG_SE_edge_type64 *y)
{
    z->a = y->a;
    z->b = y->b;
}
void LG_SE_edge2nd32_bool
    (LG_SE_edge_type32 *z, const bool *x, const LG_SE_edge_type32 *y)
{
    z->a = y->a;
    z->b = y->b;
}
void LG_SE_edge2nd64_edge
    (LG_SE_edge_type64 *z, const LG_SE_edge_type64 *x, const LG_SE_edge_type64 *y)
{
    z->a = y->a;
    z->b = y->b;
}
void LG_SE_edge2nd32_edge
    (LG_SE_edge_type32 *z, const LG_SE_edge_type32 *x, const LG_SE_edge_type32 *y)
{
    z->a = y->a;
    z->b = y->b;
}
#define EDGE2ND32_BOOL                                                           \
"void LG_SE_edge2nd32_bool                                                    \n"\
"(LG_SE_edge_type32 *z, const bool *x, const LG_SE_edge_type32 *y)            \n"\
"{                                                                            \n"\
"    //if(y->a == 0 && y->b == 0) return;                                     \n"\
"    z->a = y->a;                                                             \n"\
"    z->b = y->b;                                                             \n"\
"}"
#define EDGE2ND64_BOOL                                                           \
"void LG_SE_edge2nd64_bool                                                    \n"\
"(LG_SE_edge_type64 *z, const bool *x, const LG_SE_edge_type64 *y)            \n"\
"{                                                                            \n"\
"    //if(y->a == 0 && y->b == 0) return;                                     \n"\
"    z->a = y->a;                                                             \n"\
"    z->b = y->b;                                                             \n"\
"}"
#define EDGE2ND32_EDGE                                                           \
"void LG_SE_edge2nd32_edge                                                    \n"\
"(LG_SE_edge_type32 *z, const LG_SE_edge_type32 *x, const LG_SE_edge_type32 *y)\n"\
"{                                                                            \n"\
"    //if(y->a == 0 && y->b == 0) return;                                     \n"\
"    z->a = y->a;                                                             \n"\
"    z->b = y->b;                                                             \n"\
"}"
#define EDGE2ND64_EDGE                                                           \
"void LG_SE_edge2nd64_edge                                                    \n"\
"(LG_SE_edge_type64 *z, const LG_SE_edge_type64 *x, const LG_SE_edge_type64 *y)\n"\
"{                                                                            \n"\
"    //if(y->a == 0 && y->b == 0) return;                                     \n"\
"    z->a = y->a;                                                             \n"\
"    z->b = y->b;                                                             \n"\
"}"

int LAGr_SwapEdges
(
    // output
    LAGraph_Graph *G_new,   // A new graph with the same degree for each node
    uint64_t *pSwaps,       // Actual number of Swaps proformed
    // input: not modified
    const LAGraph_Graph G,  // Graph to be randomized.
    double loopTry,         // Percent of edges to involve per loop [0,1]
    double loopMin,         // Minimum Swaps percent per loop [0,1)
    uint64_t totSwaps,      // Desired Swaps
    uint64_t seed,          // Random Seed 
    char *msg
)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL, A_new = NULL; // n x n Adjacency Matrix 
    GrB_Vector Ai = NULL, Aj = NULL;
    // e x 1 vector, each entry is an edge.
    GrB_Vector E_vec = NULL, E_temp = NULL; 

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
    void *indices = NULL;

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
    GrB_Type M_type = NULL, E_type = NULL, Ai_type = NULL, Aj_type = NULL;
    int M_hand = 0, E_hand = 0;

    int16_t *dup_swaps = NULL;
    GrB_Vector dup_swaps_v = NULL;
    // BOOL swaps * 2 vector that holds false if an edge in the swap did not "work"

    GrB_Vector sort_h = NULL;
    GrB_Vector r_60 = NULL;
    // Only Used for coverage typecasting
    GrB_Vector temp_i = NULL;
    // count swaps 
    GrB_Index num_swaps = 0, num_attempts = 0;

    // Constants ---------------------------------------------------------------
    GrB_Vector x = NULL;

    GrB_Scalar one8 = NULL;

    GrB_Index ind_size = 0;
    
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED,
        LAGRAPH_INVALID_GRAPH, 
        "G must be undirected"
    ) ;
    // char type[LAGRAPH_MAX_NAME_LEN];
    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, 
        "G->nself_edges must be zero") ;
    LG_ASSERT (G_new != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (pSwaps != NULL, GrB_NULL_POINTER) ;
    *G_new = NULL ;

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;  

    //--------------------------------------------------------------------------
    // Extract edges from A
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new (&A_tril, GrB_BOOL, n, n)) ;
    GRB_TRY (GrB_Vector_new(&Ai, GrB_BOOL, 0)) ;
    GRB_TRY (GrB_Vector_new(&Aj, GrB_BOOL, 0)) ;
    
    // Extract lower triangular edges
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;
    GRB_TRY (GrB_Matrix_nvals(&e, A_tril)) ;
    GRB_TRY (GxB_Matrix_extractTuples_Vector(Ai, Aj, NULL, A_tril, NULL)) ;
    #ifdef COVERAGE
    if(n > 100)
    {
        // Make Ai and Aj 64 bit
        GRB_TRY (GrB_Vector_new(&temp_i, GrB_INT64, e)) ;
        GRB_TRY (GrB_assign(temp_i, NULL, NULL, Ai, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_free(&Ai)) ;
        Ai = temp_i;
        temp_i = NULL;
        GRB_TRY (GrB_Vector_new(&temp_i, GrB_INT64, e)) ;
        GRB_TRY (GrB_assign(temp_i, NULL, NULL, Aj, GrB_ALL, 0, NULL)) ;
        GRB_TRY (GrB_free(&Aj)) ;
        Aj = temp_i;
        temp_i = NULL;
    }
    #endif
    int codei = 0, codej = 0;
    GRB_TRY (GxB_Vector_type(&Ai_type, Ai)) ;
    GrB_get(Ai, &codei, GrB_EL_TYPE_CODE);
    GrB_get(Aj, &codej, GrB_EL_TYPE_CODE);

    LG_ASSERT_MSG (
        codei == codej, GrB_INVALID_VALUE,
        "extractTuples_Vector returned different types for Ai and Aj"
    ) ;
    //--------------------------------------------------------------------------
    // Initialize all operators and types
    //--------------------------------------------------------------------------
    if(codei == GrB_UINT32_CODE) // Use uint32 if possible 
    {
        GRB_TRY (GxB_Type_new(
            &lg_edge, sizeof(LG_SE_edge_type32), "LG_SE_edge_type32", EDGE_TYPE32)) ;
        GRB_TRY (GxB_Type_new(
            &lg_swap, sizeof(LG_SE_swap_type32), "LG_SE_swap_type32", SWAP_TYPE32)) ;
        GRB_TRY(GxB_BinaryOp_new(
            &hash_seed_e, (GxB_binary_function) (&LG_SE_hash_edge32),
            GrB_UINT64, lg_edge, GrB_UINT64, "LG_SE_hash_edge32", HASH_EDGE32
        )) ;
        GRB_TRY (GxB_IndexUnaryOp_new (
            &swap_pair, (GxB_index_unary_function) (&LG_SE_swap_bc32),
            lg_swap, lg_swap, GrB_BOOL, "LG_SE_swap_bc32", SWAP_BC32
        )) ;
        GRB_TRY(GxB_BinaryOp_new(
            &second_edge, (GxB_binary_function) (&LG_SE_edge2nd32_edge), 
            lg_edge, lg_edge, lg_edge, "LG_SE_edge2nd32_edge", EDGE2ND32_EDGE
        )) ;
        GRB_TRY(GxB_BinaryOp_new(
            &second_bool_edge, (GxB_binary_function) (&LG_SE_edge2nd32_bool), 
            lg_edge, GrB_BOOL, lg_edge, "LG_SE_edge2nd32_bool", EDGE2ND32_BOOL
        )) ;
    }
    else //uint64 types
    {
        GRB_TRY (GxB_Type_new(
            &lg_edge, sizeof(LG_SE_edge_type64), "LG_SE_edge_type64", EDGE_TYPE64)) ;
        GRB_TRY (GxB_Type_new(
            &lg_swap, sizeof(LG_SE_swap_type64), "LG_SE_swap_type64", SWAP_TYPE64)) ;
        GRB_TRY(GxB_BinaryOp_new(
            &hash_seed_e, (GxB_binary_function) (&LG_SE_hash_edge64),
            GrB_UINT64, lg_edge, GrB_UINT64, "LG_SE_hash_edge64", HASH_EDGE64
        )) ;
        GRB_TRY (GxB_IndexUnaryOp_new (
            &swap_pair, (GxB_index_unary_function) (&LG_SE_swap_bc64),
            lg_swap, lg_swap, GrB_BOOL, "LG_SE_swap_bc64", SWAP_BC64
        )) ;
        GRB_TRY(GxB_BinaryOp_new(
            &second_edge, (GxB_binary_function) (&LG_SE_edge2nd64_edge), 
            lg_edge, lg_edge, lg_edge, "LG_SE_edge2nd64_edge", EDGE2ND64_EDGE
        )) ;
        GRB_TRY(GxB_BinaryOp_new(
            &second_bool_edge, (GxB_binary_function) (&LG_SE_edge2nd64_bool), 
            lg_edge, GrB_BOOL, lg_edge, "LG_SE_edge2nd64_bool", EDGE2ND64_BOOL
        )) ;
    }
    
    GRB_TRY (GxB_UnaryOp_new (
        &lg_shiftland, (GxB_unary_function) (&LG_SE_shift_and),
        GrB_UINT16, GrB_UINT16, "LG_SE_shift_and", SHIFT_AND
    )) ;
    GRB_TRY(GxB_BinaryOp_new(
        &add_term_biop, (GxB_binary_function) (&LG_SE_add_term), 
        GrB_INT8, GrB_INT8, GrB_INT8, "LG_SE_add_term", ADD_TERM
    )) ;

    GRB_TRY (GxB_Monoid_terminal_new_INT8(
        &add_term_monoid, add_term_biop, (int8_t) 0, (int8_t) 2
    )) ;

    LG_SE_edge_type64 iden_second = {0,0};
    GRB_TRY (GrB_Monoid_new_UDT(
        &second_edge_monoid, second_edge, (void *) &iden_second
    )) ;

    GRB_TRY(GrB_Semiring_new(
        &plus_term_one, add_term_monoid, GrB_ONEB_INT8
    )) ;
    GRB_TRY(GrB_Semiring_new(
        &second_second_edge, second_edge_monoid, second_bool_edge
    )) ;

    //--------------------------------------------------------------------------
    // Make E Vector
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new(&E_vec, Ai_type, 2 * e)) ;
    
    // Filling out E_vec helps assign be much quicker.
    GRB_TRY (GrB_assign(
        E_vec, NULL, NULL, 0, GrB_ALL, 0, NULL));
    GrB_Index stride[] = {(GrB_Index) 0, e * 2 - 1, (GrB_Index) 2} ;

    // Shuffle i and j into E_vec.
    GRB_TRY (GrB_Vector_assign(
        E_vec, NULL, NULL, Aj, stride, GxB_STRIDE, NULL)) ;
    stride[GxB_BEGIN] = 1;
    GRB_TRY (GrB_Vector_assign(
        E_vec, NULL, NULL, Ai, stride, GxB_STRIDE, NULL)) ;

    // Reinterpret E_vec as a vector of tuples.
    GRB_TRY (GxB_Vector_unload(
        E_vec, &indices, &E_type, &e, &ind_size, &E_hand, NULL));
    e /= 2;
    GRB_TRY (GxB_Vector_load(
        E_vec, &indices, lg_edge, e, ind_size, E_hand, NULL));
    
    // Find Hash Size 
    int shift_e = 63 - (int) floor (log2 ((double) e)) ;
    uint64_t ehash_size = (1ull << (67-shift_e)) ;
    // printf("Hash Size: %ld\n", ehash_size);

    //--------------------------------------------------------------------------
    // Initialize the rest of the vectors
    //--------------------------------------------------------------------------

    // Hash values and buckets
    GRB_TRY (GrB_Vector_new(&exists, GrB_INT8, ehash_size)) ;
    GRB_TRY (GrB_Vector_new(&hashed_edges, GrB_UINT64, e)) ;
    GRB_TRY(GrB_set (exists, GxB_BITMAP | GxB_FULL, GxB_SPARSITY_CONTROL)) ;

    // Ramp 
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, e + 1)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (ramp_v, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    GrB_Index ramp_size;

    // Constants
    GRB_TRY (GrB_Scalar_new (&one8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (one8, 1)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_BOOL, e)) ;

    // Random Vector
    GRB_TRY (GrB_Vector_new(&random_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_60, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, 1ull << (64-shift_e))) ;
    GRB_TRY(GrB_set (r_permute, GxB_BITMAP, GxB_SPARSITY_CONTROL)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        random_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;
    LG_TRY(
        LAGraph_Random_Seed(random_v, seed, msg)) ;
    
    // printf("Entering loop, Good Luck:\n") ;
    while(num_swaps < totSwaps)
    {
        GrB_Index perm_size, arr_size, junk_size;
        // r_60 has the random vector shifted by some amount.
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
            r_60, NULL, NULL, GxB_BSHIFT_UINT64, random_v, -(shift_e), NULL
        )) ;
        GRB_TRY (GrB_Vector_clear(x)) ;
        GRB_TRY (GrB_Vector_resize(x, e)) ;
        GRB_TRY (GrB_Vector_assign_BOOL(
            x, NULL, NULL, true, GrB_ALL, 0, NULL)) ;
        LG_TRY (LAGraph_FastAssign_Semiring(
            r_permute, NULL, NULL, r_60, x, ramp_v, GxB_ANY_FIRSTJ_INT64,
            NULL, msg)) ;
        
        GrB_Index edges_permed = 0;
        GRB_TRY (GrB_Vector_nvals(&edges_permed, r_permute)) ;
        GRB_TRY (GrB_Vector_new(&edge_perm, GrB_BOOL, edges_permed)) ;
        GRB_TRY (GrB_Vector_extractTuples(NULL, edge_perm, r_permute, NULL)) ; 
        n_keep = LAGRAPH_MIN((int)(e * loopTry * 0.5), edges_permed / 2) ;

        // Chose only the edges we need from  edge_perm
        GRB_TRY (GrB_Vector_resize(edge_perm, n_keep * 2)) ;
        
        // Get the desired edges from the E_vec array.
        GRB_TRY (GrB_Vector_new(&M, lg_edge, n_keep * 2)) ;
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
            &new_hashed_edges, GrB_UINT64, n_keep * 2)) ;

        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
            new_hashed_edges, NULL, NULL, hash_seed_e, M, 
            ehash_size - 1ll, NULL
        )) ;
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
            hashed_edges, NULL, NULL, hash_seed_e, E_vec, 
            ehash_size - 1ll, NULL
        )) ;

        //----------------------------------------------------------------------
        // Build Hash Buckets
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Vector_new(&dup_swaps_v, GrB_INT8, n_keep * 2)) ;
        GRB_TRY (GrB_set(dup_swaps_v, GxB_BITMAP, GxB_SPARSITY_CONTROL)) ;
        
        GRB_TRY (GrB_Vector_clear(x)) ;
        GRB_TRY (GrB_Vector_resize(x, e)) ;
        GRB_TRY (GrB_Vector_assign_BOOL(
            x, NULL, NULL, true, GrB_ALL, 0, NULL)) ;
        // place a one in any bucket that coresponds to an edge currently in E
        LG_TRY (LAGraph_FastAssign_Semiring(
            exists, NULL, NULL, hashed_edges, x, ramp_v, GxB_ANY_PAIR_UINT8,
            NULL, msg
        )) ;
        
        GRB_TRY (GrB_Vector_clear(x)) ;
        GRB_TRY (GrB_Vector_resize(x, n_keep * 2)) ;
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
        // exist has to be full at this point
        

        // "Count" all of the edges that fit into each bucket. Stop counting at 
        // 2 since we will have to throw that whole bucket away anyway.
        LG_TRY (LAGraph_FastAssign_Semiring(
            exists, NULL, add_term_biop, new_hashed_edges, x, ramp_v, 
            plus_term_one, NULL, msg
        )) ;
        GRB_TRY(GrB_set (exists, GxB_BITMAP | GxB_FULL, GxB_SPARSITY_CONTROL)) ;
        // Select buckets with only one corresponding value
        GRB_TRY (GrB_Vector_select_INT8(
            exists, NULL, NULL, GrB_VALUEEQ_UINT8, exists, (int8_t) 1,
            NULL
        )) ;

        // Find each hashed edge's bucket, dup_swaps_v is 1 if exists[edge] = 1
        LG_TRY (LAGraph_FastAssign_Semiring(
            dup_swaps_v, NULL, NULL, new_hashed_edges, exists, ramp_v, 
            GxB_ANY_PAIR_INT8, GrB_DESC_T0, msg
        )) ;
        // GRB_TRY (GxB_Vector_extract_Vector(
        //     dup_swaps_v, NULL, NULL, exists, new_hashed_edges, NULL)) ;

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
        // n_keep = sum (dup_swaps_v), to count the # of 1's that now appear in
        // dup_swaps_v, which will become nvals(M) after loading M from the
        // container con.
        GRB_TRY (GrB_reduce (&n_keep, NULL, GrB_PLUS_MONOID_UINT64,
            dup_swaps_v, NULL)) ;
        con->nvals = n_keep;
        con->format = GxB_BITMAP;
        dup_swaps_v = NULL;
        GRB_TRY (GxB_load_Vector_from_Container(M, con, NULL)) ;
        // GRB_TRY (GxB_print (M, 1)) ; // for debugging; must be OK here
        GRB_TRY (GrB_free(&con)) ;
        GRB_TRY (LAGraph_FastAssign_Semiring(
            E_vec, NULL, second_edge, edge_perm, M, ramp_v, 
            second_second_edge, NULL, msg)) ;

        n_keep /= 2;

        FREE_LOOP ; // Free Matricies that have to be rebuilt

        num_swaps += n_keep ;
        LG_TRY (LAGraph_Random_Next(random_v, msg)) ;
        // printf("Made %ld swaps this loop. "
        //         "[%.3f%% of Planned, %.3f%% of edges swapped]\n"
        //         "Completed %ld swaps total. [%.3f%% of Planned]\n", 
        //      n_keep, n_keep * 100.0 / totSwaps, n_keep * 200.0 / e, num_swaps, 
        //      num_swaps * 100.0 / totSwaps) ;
        if(n_keep < (int) (loopMin * e / 2) + 1)
        {
            // printf("Too Few Swaps occured! Exiting.\n");
            break;
        }
    }
    GRB_TRY (GxB_Vector_unload(
        E_vec, (void **) &indices, &lg_edge, &e, &ind_size, &E_hand, NULL));
    GRB_TRY (GxB_Vector_load(
        E_vec, (void **) &indices, E_type, e * 2, ind_size, E_hand, NULL));
    GRB_TRY (GrB_Vector_extract(
        Aj, NULL, NULL, E_vec, stride, GxB_STRIDE, NULL));
    stride[GxB_BEGIN] = 0;
    GRB_TRY (GrB_Vector_extract(
        Ai, NULL, NULL, E_vec, stride, GxB_STRIDE, NULL));
    // Build Output Matrix
    GRB_TRY (GrB_Matrix_new(&A_new, GrB_BOOL, n, n)) ;
    GRB_TRY (GxB_Matrix_build_Scalar_Vector(A_new, Ai, Aj, one8, NULL)) ;
    GRB_TRY (GrB_eWiseAdd(
        A_new, NULL, NULL, GrB_LOR_MONOID_BOOL, A_new,A_new, GrB_DESC_T0
    )) ;
    LAGRAPH_TRY (LAGraph_New (
        G_new, &A_new, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    LG_FREE_WORK ;
    (*pSwaps) = num_swaps ;
    return (num_swaps >= totSwaps)? GrB_SUCCESS :  LAGRAPH_INSUFFICIENT_SWAPS ;
    #else
    // printf("LAGr_SwapEdges Needs GB v10\n") ;
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}
