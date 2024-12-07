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
    GrB_free (&sort_r) ;                             \
    GrB_free (&set_v) ;                             \
    GrB_free (&assign_s) ;                             \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    // start GraphBLAS and LAGraph
    bool burble = true ;               // set true for diagnostic outputs
    demo_init (burble) ;
    GrB_Matrix P = NULL;
    GrB_Vector rand_v = NULL, sort_r = NULL, set_v = NULL, assign_s = NULL,
        ramp_v = NULL;
    GrB_Index *rand_a = NULL, *ramp = NULL;
    bool *set_a = NULL;
    GrB_Index r_size = 0, ramp_size = 0, junk_size = 0;
    bool iso = false;
    LG_TRY (LAGraph_Random_Init (msg)) ;
    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).
    bool *val_of_P = NULL;
    double t = LAGraph_WallClockTime ( ) ;
    GrB_Index size = (argc > 1) ? atoll(argv [1]) : 1000 ;
    int shift_e = __builtin_clzl(size);
    GrB_Index size_p2 = (1ull << (64-shift_e));
    GrB_Index bit_mask = size_p2 - 1;
    GRB_TRY (GrB_Vector_new(&rand_v, GrB_UINT64, size)) ;
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, size + 1)) ;
    GRB_TRY (GrB_Vector_new(&sort_r, GrB_UINT64, size)) ;
    GRB_TRY (GrB_Vector_new(&set_v, GrB_BOOL, size_p2)) ;
    GRB_TRY (GrB_Vector_new(&assign_s, GrB_BOOL, size_p2)) ;
    GRB_TRY (GrB_Matrix_new(&P, GrB_BOOL, size_p2, size)) ;


    LG_TRY (LAGraph_Malloc ((void**)(&val_of_P), 1, sizeof(bool), msg)) ;
    val_of_P[0] = 1;

    GRB_TRY (GrB_Vector_assign_UINT64(
        rand_v, NULL, NULL, 0ull, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_assign_UINT64(
        ramp_v, NULL, NULL, 0ull, GrB_ALL, 0, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64(
        ramp_v, NULL, NULL, GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    GRB_TRY (GxB_Vector_unpack_Full(
        ramp_v, (void **)&ramp, &ramp_size, &iso, NULL
    )) ;
    // GRB_TRY (GrB_Vector_assign_BOOL(
    //     assign_s, NULL, NULL, 0, GrB_ALL, 0, NULL)) ;
    GRB_TRY(GrB_set (assign_s, GxB_BITMAP, GxB_SPARSITY_CONTROL) ;)
    LG_TRY (LAGraph_Random_Seed(rand_v, 1548945616ul, msg)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp1st_UINT64(
        rand_v, NULL, NULL, GrB_BAND_UINT64, bit_mask, rand_v, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    // LG_TRY (LAGraph_Vector_Print (rand_v, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // try the LAGraph_HelloWorld "algorithm"
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GxB_Vector_sort (sort_r, NULL, GrB_LT_UINT64, rand_v, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for GrB_Sort: %g sec\n", t) ;
    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GxB_Vector_unpack_Full (
        rand_v, (void **)&rand_a, &r_size, &iso, NULL
    )) ;
    LAGraph_Calloc((void **)&set_a, size_p2, sizeof(bool), msg);

    int nthreads, nthreads_outer, nthreads_inner ;
    LG_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
    nthreads = nthreads_outer * nthreads_inner ;
    printf("%d", nthreads);
    // #pragma omp parallel for num_threads(nthreads) schedule(static)
    for(int64_t i = 0; i < size; ++i)
    {
        set_a[rand_a[i]] = (bool) 1;
    }
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for Single Thread Unpack: %g sec\n", t) ;

    GRB_TRY (GxB_Vector_pack_Full (
        set_v, (void **)&set_a, 1ull << (64-shift_e), false, NULL
    )) ;
    



    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GrB_Vector_assign_BOOL(
        assign_s, NULL, NULL, 1, rand_a, size, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for Assign: %g sec\n", t) ;
    GRB_TRY (GrB_Vector_clear(assign_s)) ;
    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GxB_Matrix_pack_CSC(
        P, &ramp, &rand_a, (void**) &val_of_P, ramp_size,
        r_size, sizeof(bool), true, false, NULL
    ));
    GRB_TRY (GrB_Matrix_reduce_Monoid(
        assign_s, NULL, NULL, GxB_ANY_BOOL_MONOID, P, NULL
    ));
    GRB_TRY (GxB_Matrix_unpack_CSC(
        P, &ramp, &rand_a, (void**) &val_of_P, &ramp_size,
        &r_size, &junk_size, &iso, NULL, NULL
    ));
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for CSR Magic: %g sec\n", t) ;

    GRB_TRY (GxB_Vector_pack_Full (
        rand_v, (void **)&rand_a, r_size, iso, NULL
    )) ;
    GRB_TRY (GxB_Vector_pack_Full (
        ramp_v, (void **)&ramp, ramp_size, iso, NULL
    )) ;
    //--------------------------------------------------------------------------
    // check the results (make sure Y is a copy of G->A)
    //--------------------------------------------------------------------------
    bool isEq = 0;
    GRB_TRY (GrB_Vector_assign_BOOL(
        assign_s, assign_s, NULL, 0, GrB_ALL, 0, GrB_DESC_SC)) ;
    LG_TRY (LAGraph_Vector_IsEqual(&isEq, assign_s, set_v, msg));
    if(isEq)
        printf("TEST PASSED\n");
    else
        printf("TEST FAILED\n");

    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf ("\n===============================The result set vector:\n") ;
    // GRB_TRY (GxB_fprint(set_v, GxB_SHORT, stdout)) ;
    // GRB_TRY (GxB_fprint(assign_s, GxB_SHORT, stdout)) ;
    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    LG_TRY (LAGraph_Random_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
