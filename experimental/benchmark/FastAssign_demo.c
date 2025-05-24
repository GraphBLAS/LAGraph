//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/speed_hash_demo.c: a simple demo
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------


#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"
// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&rand_v) ;                             \
    GrB_free (&build_v) ;                             \
    GrB_free (&set_v) ;                             \
    GrB_free (&assign_s) ;                             \
    GrB_free (&fa_s) ;                             \
    GrB_free (&x) ;                             \
    GrB_free (&bool1) ;                             \
    GrB_free (&ramp) ;                             \
    LAGraph_Free ((void **)&rand_a, msg);       \
}

int main (int argc, char **argv)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;
    GrB_Vector rand_v = NULL, build_v = NULL, set_v = NULL, assign_s = NULL, 
        fa_s = NULL, x = NULL, ramp = NULL;
    GrB_Index *rand_a = NULL;
    GrB_Scalar bool1 = NULL;
    bool *set_a = NULL;
    GrB_Index r_size = 0, ramp_size = 0, junk_size = 0;
    bool iso = false;

    bool *val_of_P = NULL;
    double t = LAGraph_WallClockTime ( ) ;
    GrB_Index size = (argc > 1) ? atoll(argv [1]) : 1000 ;
    int shift_e = 63 - (int) floor (log2 ((double) size)) ;
    GrB_Index size_p2 = (1ull << (64-shift_e));
    GrB_Index bit_mask = size_p2 - 1;
    GRB_TRY (GrB_Vector_new(&rand_v, GrB_UINT64, size)) ;
    GRB_TRY (GrB_Vector_new(&x, GrB_BOOL, size)) ;
    GRB_TRY (GrB_Vector_new(&build_v, GrB_BOOL, size_p2)) ;
    GRB_TRY (GrB_Vector_new(&set_v, GrB_BOOL, size_p2)) ;
    GRB_TRY (GrB_Vector_new(&assign_s, GrB_BOOL, size_p2)) ;
    GRB_TRY (GrB_Vector_new(&fa_s, GrB_BOOL, size_p2)) ;
    GRB_TRY (GrB_Scalar_new(&bool1, GrB_BOOL));

    GrB_Type ramp_type = (size + 1 <= INT32_MAX)? GrB_UINT32: GrB_UINT64;
    GrB_IndexUnaryOp idxnum = (size + 1 <= INT32_MAX)? 
            GrB_ROWINDEX_INT32: GrB_ROWINDEX_INT64;
    GRB_TRY (GrB_Vector_new(&ramp, ramp_type, size + 1));
    GRB_TRY (GrB_assign (ramp, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_apply (ramp, NULL, NULL, idxnum, ramp, 0, NULL)) ;

    GRB_TRY (GrB_Vector_assign_UINT64(
        rand_v, NULL, NULL, 0ull, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign_BOOL(
        x, NULL, NULL, (bool) 1, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Scalar_setElement_BOOL(bool1, (bool) 1));

    GRB_TRY(GrB_set (assign_s, GxB_BITMAP, GxB_SPARSITY_CONTROL) ;)
    GRB_TRY(GrB_set (fa_s, GxB_BITMAP, GxB_SPARSITY_CONTROL) ;)
    LG_TRY (LAGraph_Random_Seed(rand_v, 1548945616ul, msg)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp1st_UINT64(
        rand_v, NULL, NULL, GrB_BAND_UINT64, bit_mask, rand_v, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to create random vector:      %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // try Methods of building a "set"
    //--------------------------------------------------------------------------

    
    // Baseline: Build
    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GxB_Vector_build_Scalar_Vector (
        build_v, rand_v, bool1, NULL)) ;
    GRB_TRY (GrB_wait(build_v, GrB_MATERIALIZE)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for Build: %g sec\n", t) ;
    t = LAGraph_WallClockTime ( ) ;

    #if 0
    GRB_TRY (GxB_Vector_unpack_Full (
        rand_v, (void **)&rand_a, &r_size, &iso, NULL
    )) ;
    // Baseline: Single Threaded random access insert
    LAGraph_Calloc((void **)&set_a, size_p2, sizeof(bool), msg);
    for(int64_t i = 0; i < size; ++i)
    {
        set_a[rand_a[i]] = (bool) 1;
    }
    t = LAGraph_WallClockTime ( ) - t ;
    GRB_TRY (GxB_Vector_pack_Full (
        set_v, (void **)&set_a, 1ull << (64-shift_e), false, NULL
    )) ;
     GRB_TRY (GxB_Vector_pack_Full (
        rand_v, (void **)&rand_a, r_size, iso, NULL
    )) ;
    printf ("Time for Single Thread Unpack: %g sec\n", t) ;
    #endif



    // Baseline: GrB_assign
    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GxB_Vector_assign_Scalar_Vector(
        assign_s, NULL, NULL, bool1, rand_v, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for GraphBLAS Assign: %g sec\n", t) ;
    
    #if GxB_IMPLEMENTATION < GxB_VERSION (10,0,0) 
    printf ("GraphBLAS version too low to test LAGraph_FastAssign\n") ;
    #else
    // FastAssign!
    t = LAGraph_WallClockTime ( ) ;
    LG_TRY (LAGraph_FastAssign_Semiring(
        fa_s, NULL, NULL, rand_v, x, ramp, GxB_ANY_PAIR_BOOL, NULL, msg
    ));
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_FastAssign: %g sec\n", t) ;
    #endif

    //--------------------------------------------------------------------------
    // check the results (Make sure that assign == build == FastAssign )
    //--------------------------------------------------------------------------
    bool flag = true, isEq = 0;
    LG_TRY (LAGraph_Vector_IsEqual(&isEq, assign_s, fa_s, msg));
    flag &= isEq;
    LG_TRY (LAGraph_Vector_IsEqual(&isEq, build_v, fa_s, msg));
    flag &= isEq;
    if(flag)
        printf("TEST PASSED\n");
    else
        printf("TEST FAILED\n");

    //--------------------------------------------------------------------------
    // print the results 
    //--------------------------------------------------------------------------

    // printf ("\n===============================The result set vector:\n") ;
    // GRB_TRY (GxB_fprint(set_v, GxB_SHORT, stdout)) ;
    // GRB_TRY (GxB_fprint(assign_s, GxB_SHORT, stdout)) ;
    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
    #else
    printf ("GraphBLAS version too low to test LAGraph_FastAssign\n") ;
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}
