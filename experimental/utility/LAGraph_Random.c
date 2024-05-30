//------------------------------------------------------------------------------
// LAGraph_Random: generate a random vector (of any sparsity structure)
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// A very simple thread-safe parallel pseudo-random nuumber generator.

// FUTURE: add LAGraph_Random_Init to LAGraph_Init,
// and added LAGraph_Random_Finalize to LAGraph_Finalize.

// LAGRAPH_V11_GENERATOR: if 1, use the simple generator in LAGraph v1.1.
// Otherwise, use xorshift64, initialized with splitmix64.
#define LAGRAPH_V11_GENERATOR 0

#include "LG_internal.h"
#include "LAGraphX.h"

void dump_state (GrB_Vector State) ;
void dump_state (GrB_Vector State)
{
    GxB_print (State, 5) ;
    uint64_t nvals ;
    GrB_Vector_nvals (&nvals, State) ;
    if (nvals > 0) 
    {
        uint64_t maxstate = 0, minstate = 0 ;
        GrB_reduce (&maxstate, NULL, GrB_MAX_MONOID_UINT64, State, NULL) ;
        GrB_reduce (&minstate, NULL, GrB_MIN_MONOID_UINT64, State, NULL) ;
        printf ("    nvals: %g min: %lu (log2 %g), max %lu (log2 %g)\n",
            (double) nvals,
            minstate, log2 ((double) minstate),
            maxstate, log2 ((double) maxstate)) ;
    }
}

//------------------------------------------------------------------------------
// LG_RAND macros
//------------------------------------------------------------------------------

// Generate the next seed
#define LG_RAND_NEXT(seed) ((seed) * 1103515245 + 12345)

//------------------------------------------------------------------------------
// global operator
//------------------------------------------------------------------------------

// These operators can be shared by all threads in a user application, and thus
// are safely declared as global objects.

GrB_UnaryOp LG_rand_next_op = NULL ;
GrB_UnaryOp LG_rand_init_op = NULL ;

//------------------------------------------------------------------------------
// LG_rand_next_f :  unary operator to construct the next state
// LG_rand_next_f2:  unary operator to construct the next state
//------------------------------------------------------------------------------

// z = f(x), where x is the old state and z is the new state.

// using xorshift, from https://en.wikipedia.org/wiki/Xorshift
// with a state of uint64_t, or xorshift64star.

// Reference: Marsaglia, George (July 2003). "Xorshift RNGs". Journal of
// Statistical Software. 8 (14).  https://doi.org/10.18637/jss.v008.i14 .

// For this random number generator, the output random number is the same
// as the state.  The default seed is given below:
#define LG_RAND_MARSAGLIA_SEED 88172645463325252LL

#if 0

    struct xorshift64_state {
        uint64_t a;
    };

    uint64_t xorshift64(struct xorshift64_state *state)
    {
            uint64_t x = state->a;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            return state->a = x;
    }

#endif

void LG_rand_next_f2 (void *z, const void *x)
{
    uint64_t state = (*((uint64_t *) x)) ;
    state ^= state << 13 ;
    state ^= state >> 7 ;
    state ^= state << 17 ;
    (*((uint64_t *) z)) = state ;
}

#define LG_RAND_NEXT_F2_DEFN                            \
"void LG_rand_next_f2 (void *z, const void *x)      \n" \
"{                                                  \n" \
"    uint64_t state = (*((uint64_t *) x)) ;         \n" \
"    state ^= state << 13 ;                         \n" \
"    state ^= state >> 7 ;                          \n" \
"    state ^= state << 17 ;                         \n" \
"    (*((uint64_t *) z)) = state ;                  \n" \
"}"

void LG_rand_next_f (void *z, const void *x)
{
    uint64_t seed = (*((uint64_t *) x)) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    (*((uint64_t *) z)) = seed ;
}

#define LG_RAND_NEXT_F_DEFN                         \
"void LG_rand_next_f (void *z, const void *x)   \n" \
"{                                              \n" \
"    uint64_t seed = (*((uint64_t *) x)) ;      \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    (*((uint64_t *) z)) = seed ;               \n" \
"}"

// Initialize the xorshift64 generator with splitmix64

// References:
//
// David Blackman and Sebastiano Vigna. Scrambled linear pseudorandom number
// generators. ACM Trans. Math. Softw., 47:1−32, 2021.
//
// Guy Steele and Sebastiano Vigna. 2021. Computationally Easy, Spectrally Good
// Multipliers for Congruential Pseudorandom Number Generators. 22 Jan.
// 2021. 23 pages. https://arxiv.org/abs/2001.05304 Revised version to appear
// in Software: Practice and Experience.  https://doi.org/10.1002/spe.3030

#if 0

    struct splitmix64_state {
        uint64_t s;
    };

    uint64_t splitmix64(struct splitmix64_state *state)
    {
        uint64_t result = (state->s += 0x9E3779B97f4A7C15);
        result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
        result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
        return result ^ (result >> 31);
    }

#endif

void LG_rand_init_func (void *z, const void *x)
{
    uint64_t state = (*((uint64_t *) x)) ;
    state += 0x9E3779B97f4A7C15 ;
    state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9 ;
    state = (state ^ (state >> 27)) * 0x94D049BB133111EB ;
    state = (state ^ (state >> 31)) ;
    if (state == 0) state = LG_RAND_MARSAGLIA_SEED ;
    (*((uint64_t *) z)) = state ;
}

#define LG_RAND_INIT_F_DEFN                                         \
"void LG_rand_init_func (void *z, const void *x)                \n" \
"{                                                              \n" \
"    uint64_t state = (*((uint64_t *) x)) ;                     \n" \
"    state += 0x9E3779B97f4A7C15 ;                              \n" \
"    state = (state ^ (state >> 30)) * 0xBF58476D1CE4E5B9 ;     \n" \
"    state = (state ^ (state >> 27)) * 0x94D049BB133111EB ;     \n" \
"    state = (state ^ (state >> 31)) ;                          \n" \
"    #define LG_RAND_MARSAGLIA_SEED 88172645463325252LL         \n" \
"    if (state == 0) state = LG_RAND_MARSAGLIA_SEED ;           \n" \
"    (*((uint64_t *) z)) = state ;                              \n" \
"}"

//------------------------------------------------------------------------------
// LAGraph_Random_Init:  create the random state operator
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_UnaryOp_free (&LG_rand_next_op) ;                   \
    GrB_UnaryOp_free (&LG_rand_init_op) ;                   \
}

int LAGraph_Random_Init (char *msg)
{
    LG_CLEAR_MSG ;
    LG_rand_next_op = NULL ;
    LG_rand_init_op = NULL ;
    #if LAGRAPH_SUITESPARSE

    #if LAGRAPH_V11_GENERATOR
    // using the generator from LAGraph v1.1
    GRB_TRY (GxB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f,
        GrB_UINT64, GrB_UINT64, "LG_rand_next_f", LG_RAND_NEXT_F_DEFN)) ;
    #else
    // using the xorshift generator from LAGraph v1.2
    GRB_TRY (GxB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f2,
        GrB_UINT64, GrB_UINT64, "LG_rand_next_f2", LG_RAND_NEXT_F2_DEFN)) ;
    #endif

    GRB_TRY (GxB_UnaryOp_new (&LG_rand_init_op, LG_rand_init_func,
        GrB_UINT64, GrB_UINT64, "LG_rand_init_func", LG_RAND_INIT_F_DEFN)) ;
    #else

    #if LAGRAPH_V11_GENERATOR
    // using the generator from LAGraph v1.1
    GRB_TRY (GrB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f,
        GrB_UINT64, GrB_UINT64)) ;
    #else
    // using the xorshift generator from LAGraph v1.2
    GRB_TRY (GrB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f2,
        GrB_UINT64, GrB_UINT64)) ;
    #endif

    GRB_TRY (GrB_UnaryOp_new (&LG_rand_init_op, LG_rand_init_func,
        GrB_UINT64, GrB_UINT64)) ;
    #endif

//  GxB_print (LG_rand_next_op, 5) ;
//  GxB_print (LG_rand_init_op, 5) ;

    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Finalize:  free the random state operator
//------------------------------------------------------------------------------

int LAGraph_Random_Finalize (char *msg)
{
    LG_CLEAR_MSG ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Seed:  create a vector of random states
//------------------------------------------------------------------------------

// Initializes a vector with random state values.  The State vector must be
// allocated on input, and should be of type GrB_UINT64.  Its sparsity
// structure is unchanged.

#undef  LG_FREE_WORK
#define LG_FREE_WORK GrB_free (&T) ;

#if defined ( COVERAGE )
// for testing only
bool random_hack = false ;
#endif

// FIXME: should this method allow the user to pass in the init and next ops?

// FIXME: rename this method:
int LAGraph_Random_Seed // construct a random state vector
(
    // input/output
    GrB_Vector State,   // GrB_UINT64 vector of random number states
    // input
    uint64_t seed,      // scalar input seed
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    GrB_Vector T = NULL ;
    LG_ASSERT (State != NULL, GrB_NULL_POINTER) ;

    // T = 1:n-1 but only for entries present in the State vector.  This
    // requires a typecast from int64 to uint64.
    GrB_Index n ;
    GRB_TRY (GrB_Vector_size (&n, State)) ;
    GRB_TRY (GrB_Vector_new (&T, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (T, NULL, NULL,
        GrB_ROWINDEX_INT64, State, 1, NULL)) ;

    #if LAGRAPH_V11_GENERATOR
    // LAGraph v1.1:
    // State = T * INT32_MAX
    GRB_TRY (GrB_apply (State, NULL, NULL, GrB_TIMES_UINT64, T,
        (uint64_t) INT32_MAX, NULL)) ;
    // State = State + seed
    GRB_TRY (GrB_apply (State, NULL, NULL, GrB_PLUS_UINT64, State, seed,
        NULL)) ;
    // State = next (State)
    GRB_TRY (GrB_Vector_apply (State, NULL, NULL, LG_rand_next_op, State,
        NULL)) ;
    #else
    // LAGraph v1.2:
    // State = next (T)
    GRB_TRY (GrB_Vector_apply (State, NULL, NULL, LG_rand_next_op, T, NULL)) ;
    // State = State + seed
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64 (State, NULL, NULL,
        GrB_PLUS_UINT64, State, seed, NULL)) ;
    // State = init (State)
    GRB_TRY (GrB_Vector_apply (State, NULL, NULL, LG_rand_init_op, State,
        NULL)) ;
    #endif

    #if defined ( COVERAGE )
    if (random_hack)
    {
        // Set all State values to 1, to break the random seed vector.
        // This is just for testing, to test algorithms that need to handle
        // extreme cases when the random number generator is non-random.
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64 (State, NULL, NULL,
            GrB_ONEB_UINT64, State, 0, NULL)) ;
    }
    #endif

//  printf ("\nseed: %" PRIu64 "\n", seed) ;
//  dump_state (State) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Next: return next vector of random seeds
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK ;

int LAGraph_Random_Next     // advance to next random vector
(
    // input/output
    GrB_Vector State,       // the sparsity pattern of State is preserved
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    LG_ASSERT (State != NULL, GrB_NULL_POINTER) ;
    // State = next (State)
    GRB_TRY (GrB_Vector_apply (State, NULL, NULL, LG_rand_next_op, State,
        NULL)) ;
//  printf ("next:\n") ;
//  dump_state (State) ;

    return (GrB_SUCCESS) ;
}

