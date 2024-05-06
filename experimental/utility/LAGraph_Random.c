//------------------------------------------------------------------------------
// LAGraph_Random: generate a random vector (of any sparsity structure)
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

// A very simple thread-safe parallel pseudo-random nuumber generator.

// FUTURE: add LAGraph_Random_Init to LAGraph_Init,
// and added LAGraph_Random_Finalize to LAGraph_Finalize.

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// LG_RAND macros
//------------------------------------------------------------------------------

// Generate the next seed
#define LG_RAND_NEXT(seed) ((seed) * 1103515245 + 12345)

// extract a random 15-bit value from a seed (no longer used)
// #define LG_RAND_15_MAX 32767
// #define LG_RAND_15(seed) (((seed)/65536) % (LG_RAND_15_MAX + 1))

//------------------------------------------------------------------------------
// global operator
//------------------------------------------------------------------------------

// This operator can be shared by all threads in a user application, and thus
// is safely declared as a global object.

GrB_UnaryOp LG_rand_next_op = NULL ;

//------------------------------------------------------------------------------
// LG_rand_next_f:  unary operator to construct the next seed
//------------------------------------------------------------------------------

// z = f(x), where x is the old seed and z is the new seed.

// FIXME: use xorshift, from https://en.wikipedia.org/wiki/Xorshift
// with a state of uint64_t, or xorshift64star

#if 0

============================================================
xorshift64:
============================================================

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

============================================================
xorshift64start:
============================================================

    but of course do not use a static variable for the state:

    /* xorshift64s, variant A_1(12,25,27) with multiplier M_32 from line 3 of
     * table 5 */
    uint64_t xorshift64star(void) {
        /* initial seed must be nonzero, don't use a static variable for the
         * state if multithreaded */
        static uint64_t x = 1;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        return x * 0x2545F4914F6CDD1DULL;
    }

#endif

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

//------------------------------------------------------------------------------
// LAGraph_Random_Init:  create the random seed operator
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                        \
{                                                           \
    GrB_UnaryOp_free (&LG_rand_next_op) ;                   \
}

int LAGraph_Random_Init (char *msg)
{
    LG_CLEAR_MSG ;
    LG_rand_next_op = NULL ;
    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f,
        GrB_UINT64, GrB_UINT64, "LG_rand_next_f", LG_RAND_NEXT_F_DEFN)) ;
    #else
    GRB_TRY (GrB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f,
        GrB_UINT64, GrB_UINT64)) ;
    #endif
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Finalize:  free the random seed operator
//------------------------------------------------------------------------------

int LAGraph_Random_Finalize (char *msg)
{
    LG_CLEAR_MSG ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Seed:  create a vector of random seeds
//------------------------------------------------------------------------------

// Initializes a vector with random seed values.  The Seed vector must be
// allocated on input, and should be of type GrB_UINT64.  Its sparsity
// structure is unchanged.

#undef  LG_FREE_WORK
#define LG_FREE_WORK GrB_free (&T) ;

#if defined ( COVERAGE )
// for testing only
bool random_hack = false ;
#endif

int LAGraph_Random_Seed     // construct a random seed vector
(
    // input/output
    GrB_Vector Seed,    // GrB_UINT64 vector of random number seeds
    // input
    uint64_t seed,      // scalar input seed
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    GrB_Vector T = NULL ;
    LG_ASSERT (Seed != NULL, GrB_NULL_POINTER) ;

    // T = 1:n but only for entries present in the Seed vector.  This
    // requires a typecast from int64 to uint64.
    GrB_Index n ;
    GRB_TRY (GrB_Vector_size (&n, Seed)) ;
    GRB_TRY (GrB_Vector_new (&T, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (T, NULL, NULL,
        GrB_ROWINDEX_INT64, Seed, 1, NULL)) ;

    // Seed = T * INT32_MAX
    GRB_TRY (GrB_apply (Seed, NULL, NULL, GrB_TIMES_UINT64, T,
        (uint64_t) INT32_MAX, NULL)) ;

    // Seed = Seed + seed
    GRB_TRY (GrB_apply (Seed, NULL, NULL, GrB_PLUS_UINT64, Seed, seed, NULL)) ;

    // Seed = next (Seed)
    GRB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;

    #if defined ( COVERAGE )
    if (random_hack)
    {
        // Set all Seed values to 1, to break the random seed vector.
        // This is just for testing, to test algorithms that need to handle
        // extreme cases when the random number generator is non-random.
        GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64 (Seed, NULL, NULL,
            GrB_ONEB_UINT64, Seed, 0, NULL)) ;
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
    GrB_Vector Seed,        // the sparsity pattern of Seed is preserved
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    LG_ASSERT (Seed != NULL, GrB_NULL_POINTER) ;
    // Seed = next (Seed)
    GRB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;
    return (GrB_SUCCESS) ;
}
