//------------------------------------------------------------------------------
// LAGraph_Random: generate a random vector (of any sparsity structure)
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// A very simple thread-safe parallel pseudo-random number generator.  These
// have a period of 2^(64)-1.  For longer periods, a better method is needed.

// In the future, the State vector could be based on a 256-bit user-defined
// type, created by LG_Random_Init below.  The LAGraph_Random_Seed and
// LAGraph_Random_Next API would not change.  Instead, they would select the
// operators based on the type of the State vector (GrB_UINT64 in the current
// method, or the 256-bit type in a future method).  Then we would need a
// single new method, say LAGraph_Random_Value, which extracts the random
// number from the State.  It would compute R = f(State), where R is GrB_UINT64
// and the State vector (with its 256-bit data type) is unchanged.  With the
// current method, R=State is implicit.

#include "LG_internal.h"

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

// using xorshift64, from https://en.wikipedia.org/wiki/Xorshift
// with a state of uint64_t.

// Reference: Marsaglia, George (July 2003). "Xorshift RNGs". Journal of
// Statistical Software. 8 (14).  https://doi.org/10.18637/jss.v008.i14 .

// For this random number generator, the output random number is the same
// as the state.

#if 0

    The default initial state is given below, but is unused here:
    #define LG_RAND_MARSAGLIA_SEED 88172645463325252LL

    struct xorshift64_state {
        uint64_t a;
    };

    // the state must be initialized to nonzero
    uint64_t xorshift64(struct xorshift64_state *state)
    {
            uint64_t x = state->a;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            return state->a = x;
    }

#endif

// return a random uint64_t; for internal use in LAGraph
uint64_t LG_Random64 (uint64_t *state)
{
    (*state) ^= (*state) << 13 ;
    (*state) ^= (*state) >> 7 ;
    (*state) ^= (*state) << 17 ;
    return (*state) ;
}

// return a random uint64_t; as a unary operator
void LG_rand_next_f2 (uint64_t *z, const uint64_t *x)
{
    uint64_t state = (*x) ;
    state ^= state << 13 ;
    state ^= state >> 7 ;
    state ^= state << 17 ;
    (*z) = state ;
}

#define LG_RAND_NEXT_F2_DEFN                                \
"void LG_rand_next_f2 (uint64_t *z, const uint64_t *x)  \n" \
"{                                                      \n" \
"    uint64_t state = (*x) ;                            \n" \
"    state ^= state << 13 ;                             \n" \
"    state ^= state >> 7 ;                              \n" \
"    state ^= state << 17 ;                             \n" \
"    (*z) = state ;                                     \n" \
"}"

// From these references, the recommendation is to create the initial state of
// a random number generator with an entirely different random number
// generator.  splitmix64 is recommended, so we initialize the State(i) with
// splitmix64 (i+seed).  The method cannot return a value of zero, so it is
// suitable as a seed for the xorshift64 generator, above.

// References:
//
// David Blackman and Sebastiano Vigna. Scrambled linear pseudorandom number
// generators. ACM Trans. Math. Softw., 47:1−32, 2021.
//
// Steele GL, Vigna S. Computationally easy, spectrally good multipliers for
// congruential pseudorandom number generators.  Software: Practice and
// Experience 2022; 52(2): 443–458. https://doi.org/10.1002/spe.3030
//
// Guy L. Steele, Doug Lea, and Christine H. Flood. 2014. Fast splittable
// pseudorandom number generators. SIGPLAN Not. 49, 10 (October 2014), 453–472.
// https://doi.org/10.1145/2714064.2660195
//
// The splitmix64 below method is the mix64variant13 in the above paper.

#define GOLDEN_GAMMA 0x9E3779B97F4A7C15LL

#if 0

    struct splitmix64_state {
        uint64_t s;
    };

    uint64_t splitmix64(struct splitmix64_state *state)
    {
        uint64_t result = (state->s += 0x9E3779B97F4A7C15LL) ;
        result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9LL ;
        result = (result ^ (result >> 27)) * 0x94D049BB133111EBLL ;
        return result ^ (result >> 31);
    }

#endif

// The init function computes z = splitmix64 (i + seed), but it does not
// advance the seed value on return.
void LG_rand_init_func (uint64_t *z, const void *x,
    GrB_Index i, GrB_Index j, const uint64_t *seed)
{
    uint64_t state = i + (*seed) ;
    uint64_t result = (state += GOLDEN_GAMMA) ;
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9LL ;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EBLL ;
    result = (result ^ (result >> 31)) ;
    (*z) = result ;
}

#define LG_RAND_INIT_F_DEFN                                         \
"void LG_rand_init_func (uint64_t *z, const void *x,            \n" \
"    GrB_Index i, GrB_Index j, const uint64_t *seed)            \n" \
"{                                                              \n" \
"   uint64_t state = i + (*seed) ;                              \n" \
"   uint64_t result = (state += 0x9E3779B97F4A7C15LL) ;         \n" \
"   result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9LL ; \n" \
"   result = (result ^ (result >> 27)) * 0x94D049BB133111EBLL ; \n" \
"   result = (result ^ (result >> 31)) ;                        \n" \
"   (*z) = result ;                                             \n" \
"}"

//------------------------------------------------------------------------------
// LG_Random_Init:  create the random state operator
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_UnaryOp_free (&LG_rand_next_op) ;                   \
    GrB_IndexUnaryOp_free (&LG_rand_init_op) ;              \
}

#if !LAGRAPH_SUITESPARSE
typedef void (*grb_unary_function)  (void *, const void *) ;
typedef void (*grb_index_unary_function)
(
    void *z,            // output value z, of type ztype
    const void *x,      // input value x of type xtype; value of v(i) or A(i,j)
    GrB_Index i,        // row index of A(i,j)
    GrB_Index j,        // column index of A(i,j), or zero for v(i)
    const void *y       // input scalar y
) ;
#endif

int LG_Random_Init (char *msg)
{
    LG_CLEAR_MSG ;
    LG_FREE_WORK ; // free the two ops in case LG_Random_Init is called twice
    LG_rand_next_op = NULL ;
    LG_rand_init_op = NULL ;

    #if LAGRAPH_SUITESPARSE
    {
        // give SuiteSparse:GraphBLAS the strings that define the functions
        GRB_TRY (GxB_UnaryOp_new (&LG_rand_next_op,
            (GxB_unary_function) LG_rand_next_f2,
            GrB_UINT64, GrB_UINT64,
            "LG_rand_next_f2", LG_RAND_NEXT_F2_DEFN)) ;
        GRB_TRY (GxB_IndexUnaryOp_new (&LG_rand_init_op,
            (GxB_index_unary_function) LG_rand_init_func,
            GrB_UINT64, GrB_UINT64, GrB_UINT64,
            "LG_rand_init_func", LG_RAND_INIT_F_DEFN)) ;
    }
    #else
    {
        // vanilla GraphBLAS, no strings to define the new operators
        GRB_TRY (GrB_UnaryOp_new (&LG_rand_next_op,
            (grb_unary_function) LG_rand_next_f2,
            GrB_UINT64, GrB_UINT64)) ;
        GRB_TRY (GrB_IndexUnaryOp_new (&LG_rand_init_op,
            (grb_index_unary_function) LG_rand_init_func,
            GrB_UINT64, GrB_UINT64, GrB_UINT64)) ;
    }
    #endif

    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LG_Random_Finalize:  free the random state operator
//------------------------------------------------------------------------------

int LG_Random_Finalize (char *msg)
{
    LG_CLEAR_MSG ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Seed:  create a vector of random states
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK ;

// Initializes a vector with random state values.  The State vector must be
// allocated on input, and should be of type GrB_UINT64.  Its sparsity
// structure is unchanged.

#if defined ( COVERAGE )
// for testing only
bool random_hack = false ;
#endif

int LAGraph_Random_Seed // construct a random State vector
(
    // input/output:
    GrB_Vector State,   // vector of random number States, normally GrB_UINT64
    // input:
    uint64_t seed,      // scalar input seed
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    LG_ASSERT (State != NULL, GrB_NULL_POINTER) ;

    // State (i) = splitmix64 (i + seed) for all prior entries in State
    GRB_TRY (GrB_apply (State, NULL, NULL, LG_rand_init_op, State, seed,
        NULL)) ;

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

    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Next: return next vector of random seeds
//------------------------------------------------------------------------------

int LAGraph_Random_Next     // advance to next random vector
(
    // input/output:
    GrB_Vector State,   // vector of random number States, normally GrB_UINT64
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    LG_ASSERT (State != NULL, GrB_NULL_POINTER) ;
    // State = xorshift64 (State)
    GRB_TRY (GrB_apply (State, NULL, NULL, LG_rand_next_op, State, NULL)) ;
    return (GrB_SUCCESS) ;
}

