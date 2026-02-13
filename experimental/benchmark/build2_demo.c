//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/build2_demo.c: benchmark build
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

// usage:  ./build2_demo m n nvals seed
// to build a random m-by-n matrix with nvals entries
// If not present, m=4, n=6, nvals=10, seed=1 is the default.

// This main program makes use of supporting utilities in
// src/benchmark/LAGraph_demo.h and src/utility/LG_internal.h.

// GraphBLAS v10 is required for this demo.

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&Mod) ;                           \
    GrB_free (&A) ;                             \
    GrB_free (&B) ;                             \
    GrB_free (&I) ;                             \
    GrB_free (&J) ;                             \
    GrB_free (&X) ;                             \
    GrB_free (&State) ;                         \
}

//------------------------------------------------------------------------------
// mod function for uint64: z = x % y
//------------------------------------------------------------------------------

void LG_randmat_mod (void *z, const void *x, const void *y)
{
    uint64_t a = (*((uint64_t *) x)) ;
    uint64_t b = (*((uint64_t *) y)) ;
    (*((uint64_t *) z)) = a % b ;
}

#define MOD_FUNCTION_DEFN                                           \
"void LG_randmat_mod (void *z, const void *x, const void *y)    \n" \
"{                                                              \n" \
"    uint64_t a = (*((uint64_t *) x)) ;                         \n" \
"    uint64_t b = (*((uint64_t *) y)) ;                         \n" \
"    (*((uint64_t *) z)) = a % b ;                              \n" \
"}"

int main (int argc, char **argv)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    GrB_BinaryOp Mod = NULL ;
    GrB_Matrix A = NULL, B = NULL ;
    GrB_Vector I = NULL, J = NULL, X = NULL, State = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;              // set true for diagnostic outputs
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // get inputs
    //--------------------------------------------------------------------------

    GrB_Index nrows = 4, ncols = 6, nvals = 10, seed = 1 ;
    if (argc > 1) nrows = atoi (argv [1]) ;
    if (argc > 2) ncols = atoi (argv [2]) ;
    if (argc > 3) nvals = atoi (argv [3]) ;
    if (argc > 4) nvals = atoi (argv [3]) ;
    printf ("Generating a random %lu-by-%lu matrix with %lu entries"
        " (seed: %lu)\n", nrows, ncols, nvals, seed) ;

    GRB_TRY (GxB_BinaryOp_new (&Mod, LG_randmat_mod,
        GrB_UINT64, GrB_UINT64, GrB_UINT64,
        "LG_randmat_mod", MOD_FUNCTION_DEFN)) ;

    //--------------------------------------------------------------------------
    // construct the random tuples
    //--------------------------------------------------------------------------

    double t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GrB_Vector_new (&State, GrB_UINT64, nvals)) ;
    GRB_TRY (GrB_assign (State, NULL, NULL, 0, GrB_ALL, nvals, NULL)) ;
    LG_TRY (LAGraph_Random_Seed (State, seed, msg)) ;

    // I = mod (State, nrows) ;
    GRB_TRY (GrB_Vector_new (&I, GrB_UINT64, nvals)) ;
    GRB_TRY (GrB_apply (I, NULL, NULL, Mod, State, nrows, NULL)) ;

    // State = next (State)
    LG_TRY (LAGraph_Random_Next (State, msg)) ;

    // J = mod (State, ncols) ;
    GRB_TRY (GrB_Vector_new (&J, GrB_UINT64, nvals)) ;
    GRB_TRY (GrB_apply (J, NULL, NULL, Mod, State, ncols, NULL)) ;

    // State = next (State)
    LG_TRY (LAGraph_Random_Next (State, msg)) ;

    // X = (double) State
    GRB_TRY (GrB_Vector_new (&X, GrB_FP64, nvals)) ;
    GRB_TRY (GrB_assign (X, NULL, NULL, State, GrB_ALL, nvals, NULL)) ;
    GrB_free (&State) ;

    // X = X / (double) UINT64_MAX
    GRB_TRY (GrB_apply (X, NULL, NULL, GrB_DIV_FP64, X, (double) UINT64_MAX,
        NULL)) ;

    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to construct random tuples: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // build the matrix
    //--------------------------------------------------------------------------

    for (int k = 0 ; k < 3 ; k++)
    {
        t = LAGraph_WallClockTime ( ) ;
        GRB_TRY (GrB_Matrix_new (&A, GrB_FP64, nrows, ncols)) ;
        GRB_TRY (GrB_Matrix_set_INT32 (A, GxB_HYPERSPARSE,
            GxB_SPARSITY_CONTROL)) ;
        GRB_TRY (GxB_Matrix_build_Vector (A, I, J, X, GrB_PLUS_FP64, NULL)) ;
        t = LAGraph_WallClockTime ( ) - t ;
        printf ("Time for build (%d):         %g sec\n", k, t) ;

        // TODO: check the results, CPU vs GPU
        // LG_TRY (LAGraph_Matrix_IsEqual (&isequal, A, B, msg)) ;

        GRB_TRY (GxB_print (A, 1)) ;

        t = LAGraph_WallClockTime ( ) ;
        double sum = 0 ;
        GRB_TRY (GrB_Matrix_reduce_FP64 (&sum, NULL, GrB_PLUS_MONOID_FP64, A,
            NULL)) ;
        t = LAGraph_WallClockTime ( ) - t ;
        printf ("sum %g\n", sum) ;
        printf ("Time for reduce (%d)         %g sec\n", k, t) ;
        GrB_Matrix_free (&A) ;
    }

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
#endif
    return (GrB_SUCCESS) ;
}

