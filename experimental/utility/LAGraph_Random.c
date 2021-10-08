//------------------------------------------------------------------------------
// LAGraph_Random: generate a random vector (of any sparsity pattern)
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// Contributed by Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// A simple thread-safe parallel pseudo-random nuumber generator.

// FUTURE: add LAGraph_Random_Init to LAGraph_Init,
// and added LAGraph_Random_Finalize to LAGraph_Finalize.

#include "LG_internal.h"

//------------------------------------------------------------------------------
// LG_RAND macros
//------------------------------------------------------------------------------

// Generate the next seed, and extract a random 15-bit value from a seed.

#define LG_RAND_NEXT(seed) ((seed) * 1103515245 + 12345)
#define LG_RAND_15_MAX 32767 
#define LG_RAND_15(seed) (((seed)/65536) % (LG_RAND_15_MAX + 1))

//------------------------------------------------------------------------------
// global operators
//------------------------------------------------------------------------------

// These can be shared by all threads in a user application, and thus are
// safely declared as global objects.

GrB_UnaryOp LG_rand_next_op = NULL ;
GrB_UnaryOp LG_rand_iget_op = NULL ;
GrB_UnaryOp LG_rand_xget32_op = NULL ;
GrB_UnaryOp LG_rand_xget64_op = NULL ;

//------------------------------------------------------------------------------
// LG_rand_next_op:  unary operator to construct the next seed
//------------------------------------------------------------------------------

// z = f(x), where x is the old seed and z is the new seed.

void LG_rand_next_f (void *z, const void *x)
{
    uint64_t seed = (uint64_t) (*((int64_t *) x)) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    seed = LG_RAND_NEXT (seed) ;
    (*((int64_t *) z)) = (int64_t) seed ;
}

//------------------------------------------------------------------------------
// LG_rand_iget_f:  unary op to construct get a random integer from the seed
//------------------------------------------------------------------------------

// z = f(x), where x is a random seed, and z is an signed 64-bit
// pseudo-random number constructed from the seed.

void LG_rand_iget_f (void *z, const void *x)
{
    uint64_t seed = (uint64_t) (*((int64_t *) x)) ;
    uint64_t r =             LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    (*((int64_t *) z)) = (int64_t) r ;
}

//------------------------------------------------------------------------------
// LG_rand_xget64_f:  unary op to construct get a random double from the seed
//------------------------------------------------------------------------------

// z = f(x), where x is a random seed, and z is a double precision
// pseudo-random number constructed from the seed, in the range 0 to 1.

void LG_rand_xget64_f (void *z, const void *x)
{
    uint64_t seed = (uint64_t) (*((int64_t *) x)) ;
    uint64_t r =             LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    (*((double *) z)) = ((double) r) / ((double) UINT64_MAX) ; 
}

//------------------------------------------------------------------------------
// LG_rand_xget32_f:  unary op to construct get a random float from the seed
//------------------------------------------------------------------------------

// z = f(x), where x is a random seed, and z is a double precision
// pseudo-random number constructed from the seed, in the range 0 to 1.

void LG_rand_xget32_f (void *z, const void *x)
{
    uint64_t seed = (uint64_t) (*((int64_t *) x)) ;
    uint64_t r =             LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    r = LG_RAND_15_MAX * r + LG_RAND_15 (seed) ; seed = LG_RAND_NEXT (seed) ;
    (*((float *) z)) = ((float) r) / ((float) UINT64_MAX) ; 
}

//------------------------------------------------------------------------------
// LAGraph_Random_Init:  create the random seed operators
//------------------------------------------------------------------------------

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK                                   \
{                                                           \
    GrB_UnaryOp_free (&LG_rand_next_op) ;                   \
    GrB_UnaryOp_free (&LG_rand_iget_op) ;                   \
    GrB_UnaryOp_free (&LG_rand_xget32_op) ;                 \
    GrB_UnaryOp_free (&LG_rand_xget64_op) ;                 \
}

int LAGraph_Random_Init (char *msg)
{
    LG_CLEAR_MSG ;
    LG_rand_next_op = NULL ;
    LG_rand_iget_op = NULL ;
    LG_rand_xget32_op = NULL ;
    LG_rand_xget64_op = NULL ;
    GrB_TRY (GrB_UnaryOp_new (&LG_rand_next_op, LG_rand_next_f,
        GrB_INT64, GrB_INT64)) ;
    GrB_TRY (GrB_UnaryOp_new (&LG_rand_iget_op, LG_rand_iget_f,
        GrB_INT64, GrB_INT64)) ;
    GrB_TRY (GrB_UnaryOp_new (&LG_rand_xget32_op, LG_rand_xget32_f,
        GrB_FP32, GrB_INT64)) ;
    GrB_TRY (GrB_UnaryOp_new (&LG_rand_xget64_op, LG_rand_xget64_f,
        GrB_FP64, GrB_INT64)) ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Finalize:  free the random seed operators
//------------------------------------------------------------------------------

int LAGraph_Random_Finalize (char *msg)
{
    LG_CLEAR_MSG ;
    LAGraph_FREE_WORK ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Seed:  create a vector of random seeds
//------------------------------------------------------------------------------

// Initializes a vector with random seed values.  The Seed vector must be
// allocated on input.  Its sparsity pattern is unchanged.

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK ;

int LAGraph_Random_Seed     // construct a random seed vector
(
    // input/output
    GrB_Vector Seed,    // vector of random number seeds
    // input
    int64_t seed,       // scalar input seed
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    if (Seed == NULL) return (-1) ;

    #if LG_SUITESPARSE

        #if (GxB_IMPLEMENTATION >= GxB_VERSION (5,2,0))

            // Seed = (0:n-1) + seed
            // requires v2.0 C API
            GrB_TRY (GrB_Vector_apply_IndexOp_INT64 (Seed, NULL, NULL,
                GrB_ROWINDEX_INT64, Seed, seed, NULL)) ;

        #else

            // Seed <Seed> = seed
            GrB_Index n ;
            GrB_TRY (GrB_Vector_size (&n, Seed)) ;
            GrB_TRY (GrB_Vector_assign_INT64 (Seed, Seed, NULL, seed,
                GrB_ALL, n, GrB_DESC_S)) ;
            // Seed += 0:n-1
            GrB_TRY (GrB_Vector_apply (Seed, NULL, GrB_PLUS_INT64,
                GxB_POSITIONI_INT64, Seed, NULL)) ;

        #endif

    #else

            // Seed = (0:n-1) + seed
            // requires v2.0 C API
            GrB_TRY (GrB_Vector_apply_IndexOp_INT64 (Seed, NULL, NULL,
                GrB_ROWINDEX_INT64, Seed, seed, NULL)) ;

    #endif

    // Seed = next (Seed)
    GrB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_Next: return next vector of random seeds
//------------------------------------------------------------------------------

int LAGraph_Random_Next     // random int64 vector of seeds
(
    // input/output
    GrB_Vector Seed,
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    if (Seed == NULL) return (-1) ;
    // Seed = next (Seed)
    GrB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_INT64: return a vector of random int64 integers
//------------------------------------------------------------------------------

// The sparsity pattern of the result X is the same as the Seed vector.
// The Seed vector is updated to advance to the next set of seeds.

int LAGraph_Random_INT64    // random int64 vector
(
    // output
    GrB_Vector X,       // already allocated on input
    // input/output
    GrB_Vector Seed,
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    if (Seed == NULL || X == NULL) return (-1) ;
    // X = iget (Seed)
    GrB_TRY (GrB_Vector_apply (X,    NULL, NULL, LG_rand_iget_op, Seed, NULL)) ;
    // Seed = next (Seed)
    GrB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_Random_FP64: return a vector of random doubles, 0 to 1 inclusive
//------------------------------------------------------------------------------

// The sparsity pattern of the result X is the same as the Seed vector.
// The Seed vector is updated to advance to the next set of seeds.

GrB_Info LAGraph_Random_FP64    // random double vector
(
    // output
    GrB_Vector X,       // already allocated on input
    // input/output
    GrB_Vector Seed,
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    if (Seed == NULL || X == NULL) return (-1) ;
    // X = xget (Seed)
    GrB_TRY (GrB_Vector_apply (X,    NULL, NULL, LG_rand_xget64_op, Seed,
        NULL)) ;
    // Seed = next (Seed)
    GrB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;
    return (0) ;
}


//------------------------------------------------------------------------------
// LAGraph_Random_FP32: return a vector of random floats, 0 to 1 inclusive
//------------------------------------------------------------------------------

// The sparsity pattern of the result X is the same as the Seed vector.
// The Seed vector is updated to advance to the next set of seeds.

GrB_Info LAGraph_Random_FP32    // random float vector
(
    // output
    GrB_Vector X,       // already allocated on input
    // input/output
    GrB_Vector Seed,
    char *msg
)
{
    // check inputs
    LG_CLEAR_MSG ;
    if (Seed == NULL || X == NULL) return (-1) ;
    // X = xget (Seed)
    GrB_TRY (GrB_Vector_apply (X,    NULL, NULL, LG_rand_xget32_op, Seed,
        NULL)) ;
    // Seed = next (Seed)
    GrB_TRY (GrB_Vector_apply (Seed, NULL, NULL, LG_rand_next_op, Seed, NULL)) ;
    return (0) ;
}

