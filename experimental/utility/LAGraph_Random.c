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

// FIXME: ready for src?

// FIXME: add LAGraph_Random_Init to LAGraph_Init,
// and added LAGraph_Random_Finalize to LAGraph_Finalize.

// FIXME: is the new init function more complicated than it needs to be?

// LAGRAPH_V11_GENERATOR: if 1, use the simple generator in LAGraph v1.1.
// Otherwise, use xorshift64, initialized with splitmix64.
#define LAGRAPH_V11_GENERATOR 0

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// global operator
//------------------------------------------------------------------------------

// These operators can be shared by all threads in a user application, and thus
// are safely declared as global objects.

GrB_UnaryOp LG_rand_next_op = NULL ;
GrB_IndexUnaryOp LG_rand_init_op = NULL ;

//------------------------------------------------------------------------------
// unary and index-unary ops to construct the first and next states
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

#define LG_RAND_NEXT(seed) ((seed) * 1103515245 + 12345)
void LG_rand_next_f1 (void *z, const void *x)
{
    uint64_t seed = (*((uint64_t *) x)) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    (*((uint64_t *) z)) = seed ;
}

#define LG_RAND_NEXT_F1_DEFN                        \
"void LG_rand_next_f1 (void *z, const void *x)  \n" \
"{                                              \n" \
"    uint64_t seed = (*((uint64_t *) x)) ;      \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    seed = ((seed) * 1103515245 + 12345) ;     \n" \
"    (*((uint64_t *) z)) = seed ;               \n" \
"}"

// From these references, the recommendation is to create the initial state of
// a random number generator with an entirely different random number
// generator.  splitmix64 is recommended, but here we initialize the State(i)
// with xorshift (i+1) to get a good start, add the scalar seed, and then
// randomize the State with splitmix64.

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

// The init function computes z = splitmix64 (xorshift (i+1) + seed)
void LG_rand_init_func (void *z, const void *x, GrB_Index i,
    GrB_Index j, const void *y)
{
    // state = xorshift64 (i+1) + seed
    uint64_t state = i + 1 ;
    state ^= state << 13 ;
    state ^= state >> 7 ;
    state ^= state << 17 ;
    uint64_t seed = (*((uint64_t *) y)) ;
    state += seed ;
    // result = shiftmix64 (state)
    uint64_t result = (state += 0x9E3779B97f4A7C15) ;
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9 ;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB ;
    result = (result ^ (result >> 31)) ;
    // this is a precaution against the unlikely event that state is zero:
    if (result == 0) result = LG_RAND_MARSAGLIA_SEED ;
    // return the result
    (*((uint64_t *) z)) = result ;
}

#define LG_RAND_INIT_F_DEFN                                         \
"void LG_rand_init_func (void *z, const void *x, GrB_Index i,   \n" \
"    GrB_Index j, const void *y)                                \n" \
"{                                                              \n" \
"   uint64_t state = i + 1 ;                                    \n" \
"   state ^= state << 13 ;                                      \n" \
"   state ^= state >> 7 ;                                       \n" \
"   state ^= state << 17 ;                                      \n" \
"   uint64_t seed = (*((uint64_t *) y)) ;                       \n" \
"   state += seed ;                                             \n" \
"   uint64_t result = (state += 0x9E3779B97f4A7C15) ;           \n" \
"   result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9 ;   \n" \
"   result = (result ^ (result >> 27)) * 0x94D049BB133111EB ;   \n" \
"   result = (result ^ (result >> 31)) ;                        \n" \
"   if (result == 0) result = 88172645463325252LL ;             \n" \
"   (*((uint64_t *) z)) = result ;                              \n" \
"}"

//------------------------------------------------------------------------------
// LAGraph_Random_Init:  create the random state operator
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_UnaryOp_free (&LG_rand_next_op) ;                   \
    GrB_IndexUnaryOp_free (&LG_rand_init_op) ;              \
}

int LAGraph_Random_Init (char *msg)
{
    LG_CLEAR_MSG ;
    LG_rand_next_op = NULL ;
    LG_rand_init_op = NULL ;

    #if LAGRAPH_SUITESPARSE
    {
        // give SuiteSparse:GraphBLAS the strings that define the functions
        #if LAGRAPH_V11_GENERATOR
        {
            // using the generator from LAGraph v1.1
            GRB_TRY (GxB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f1,
                GrB_UINT64, GrB_UINT64,
                "LG_rand_next_f1", LG_RAND_NEXT_F1_DEFN)) ;
        }
        #else
        {
            // using the xorshift generator from LAGraph v1.2
            GRB_TRY (GxB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f2,
                GrB_UINT64, GrB_UINT64,
                "LG_rand_next_f2", LG_RAND_NEXT_F2_DEFN)) ;
            GRB_TRY (GxB_IndexUnaryOp_new (&LG_rand_init_op, LG_rand_init_func,
                GrB_UINT64, GrB_UINT64, GrB_UINT64,
                "LG_rand_init_func", LG_RAND_INIT_F_DEFN)) ;
        }
        #endif
    }
    #else
    {
        // vanilla GraphBLAS, no strings to define the new operators
        #if LAGRAPH_V11_GENERATOR
        {
            // using the generator from LAGraph v1.1
            GRB_TRY (GrB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f1,
                GrB_UINT64, GrB_UINT64)) ;
        }
        #else
        {
            // using the xorshift generator from LAGraph v1.2
            GRB_TRY (GrB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f2,
                GrB_UINT64, GrB_UINT64)) ;
            GRB_TRY (GrB_IndexUnaryOp_new (&LG_rand_init_op, LG_rand_init_func,
                GrB_UINT64, GrB_UINT64, GrB_UINT64) ;
        }
        #endif
    }
    #endif

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
    GrB_Vector T = NULL ;
    LG_CLEAR_MSG ;
    LG_ASSERT (State != NULL, GrB_NULL_POINTER) ;

    #if LAGRAPH_V11_GENERATOR
    {
        // LAGraph v1.1:
        // T = 1:n but only for the entries present in the State vector.  This
        // requires a typecast from int64 to uint64.
        GrB_Index n ;
        GRB_TRY (GrB_Vector_size (&n, State)) ;
        GRB_TRY (GrB_Vector_new (&T, GrB_UINT64, n)) ;
        GRB_TRY (GrB_apply (T, NULL, NULL,
            GrB_ROWINDEX_INT64, State, 1, NULL)) ;
        // State = T * INT32_MAX
        GRB_TRY (GrB_apply (State, NULL, NULL, GrB_TIMES_UINT64, T,
            (uint64_t) INT32_MAX, NULL)) ;
        // State = State + seed
        GRB_TRY (GrB_apply (State, NULL, NULL, GrB_PLUS_UINT64, State, seed,
            NULL)) ;
        // State = next (State)
        GRB_TRY (GrB_apply (State, NULL, NULL, LG_rand_next_op, State, NULL)) ;
    }
    #else
    {
        // LAGraph v1.2:
        // State = splitmix64 (xorshift64 (i+1) + seed)
        GRB_TRY (GrB_apply (State, NULL, NULL, LG_rand_init_op, State, seed,
            NULL)) ;
    }
    #endif

    #if defined ( COVERAGE )
    if (random_hack)
    {
        // Set all State values to 1, to break the random seed vector.
        // This is just for testing, to test algorithms that need to handle
        // extreme cases when the random number generator is non-random.
        GRB_TRY (GrB_apply (State, NULL, NULL, GrB_ONEB_UINT64, State, 0,
            NULL)) ;
    }
    #endif

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
    GRB_TRY (GrB_apply (State, NULL, NULL, LG_rand_next_op, State, NULL)) ;
    return (GrB_SUCCESS) ;
}

