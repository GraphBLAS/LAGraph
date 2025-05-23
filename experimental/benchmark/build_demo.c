//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/build_demo.c: benchmark build & setElement
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
    LAGraph_Free ((void **) &I2, NULL) ;        \
    LAGraph_Free ((void **) &J2, NULL) ;        \
    LAGraph_Free ((void **) &X2, NULL) ;        \
    GrB_free (&A) ;                             \
    GrB_free (&I) ;                             \
    GrB_free (&J) ;                             \
    GrB_free (&X) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Vector I = NULL, J = NULL, X = NULL ;
    void *I2 = NULL, *J2 = NULL, *X2 = NULL ;
    uint64_t I2_size = 0, J2_size = 0, X2_size = 0,
             I2_len  = 0, J2_len  = 0, X2_len  = 0 ;
    int I2_handling = 0, J2_handling = 0, X2_handling = 0 ;
    GrB_Type I2_type = NULL, J2_type = NULL, X2_type = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).

    double t = LAGraph_WallClockTime ( ) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        false,      // make the graph undirected, if true
        false,      // remove self-edges, if true
        false,      // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // ensure G->A is non-iso
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GrB_Matrix_set_INT32 (G->A, false, GxB_ISO)) ;
    GRB_TRY (GrB_Matrix_setElement_INT32 (G->A, 0, 0, 0)) ;
    GRB_TRY (GrB_wait (G->A, GrB_MATERIALIZE)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to make A non-iso: %g sec\n", t) ;

//  printf ("\n==========================The input graph matrix G:\n") ;
//  LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // extract all tuples from G->A
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GrB_Vector_new (&I, GrB_UINT64, 0)) ;
    GRB_TRY (GrB_Vector_new (&J, GrB_UINT64, 0)) ;
    GRB_TRY (GrB_Vector_new (&X, GrB_UINT64, 0)) ;
    GRB_TRY (GxB_Matrix_extractTuples_Vector (I, J, X, G->A, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for extractTuples: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // build a copy of G->A using GrB_build
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    GrB_Index nrows, ncols ;
    GRB_TRY (GrB_Matrix_nrows (&nrows, G->A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    GrB_Type atype ;
    GRB_TRY (GxB_Matrix_type (&atype, G->A)) ;
    GRB_TRY (GrB_Matrix_new (&A, atype, nrows, ncols)) ;
    GRB_TRY (GxB_Matrix_build_Vector (A, I, J, X, GxB_IGNORE_DUP, NULL)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for build:         %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results (make sure A is a copy of G->A)
    //--------------------------------------------------------------------------

    bool isequal ;
    t = LAGraph_WallClockTime ( ) ;
    LG_TRY (LAGraph_Matrix_IsEqual (&isequal, A, G->A, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check:          %g sec\n", t) ;
    LG_ASSERT (isequal, GrB_INVALID_VALUE) ;

    //--------------------------------------------------------------------------
    // build a copy of G->A using setElement
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Matrix_new (&A, atype, nrows, ncols)) ;
    GRB_TRY (GxB_Vector_unload (I, &I2, &I2_type, &I2_len, &I2_size,
        &I2_handling, NULL)) ;
    GRB_TRY (GxB_Vector_unload (J, &J2, &J2_type, &J2_len, &J2_size,
        &J2_handling, NULL)) ;
    GRB_TRY (GxB_Vector_unload (X, &X2, &X2_type, &X2_len, &X2_size,
        &X2_handling, NULL)) ;

    uint64_t nvals = I2_len ;
    LG_ASSERT (nvals == J2_len && nvals == X2_len, GrB_INVALID_VALUE) ;

    uint32_t *I32 = (I2_type == GrB_UINT32) ? (uint32_t *) I2 : NULL ;
    uint64_t *I64 = (I2_type == GrB_UINT32) ? NULL : (uint64_t *) I2 ;

    uint32_t *J32 = (J2_type == GrB_UINT32) ? (uint32_t *) J2 : NULL ;
    uint64_t *J64 = (J2_type == GrB_UINT32) ? NULL : (uint64_t *) J2 ;

    #define BUILD_WITH_SETEL(xtype)                                     \
    {                                                                   \
        xtype *XX = (xtype *) X2 ;                                      \
        if (I32 && J32)                                                 \
        {                                                               \
            for (int64_t k = 0 ; k < nvals ; k++)                       \
            {                                                           \
                GrB_Index i = I32 [k] ;                                 \
                GrB_Index j = J32 [k] ;                                 \
                GRB_TRY (GrB_Matrix_setElement (A, XX [k], i, j)) ;     \
            }                                                           \
        }                                                               \
        else if (I32 && J64)                                            \
        {                                                               \
            for (int64_t k = 0 ; k < nvals ; k++)                       \
            {                                                           \
                GrB_Index i = I32 [k] ;                                 \
                GrB_Index j = J64 [k] ;                                 \
                GRB_TRY (GrB_Matrix_setElement (A, XX [k], i, j)) ;     \
            }                                                           \
        }                                                               \
        else if (I64 && J32)                                            \
        {                                                               \
            for (int64_t k = 0 ; k < nvals ; k++)                       \
            {                                                           \
                GrB_Index i = I64 [k] ;                                 \
                GrB_Index j = J32 [k] ;                                 \
                GRB_TRY (GrB_Matrix_setElement (A, XX [k], i, j)) ;     \
            }                                                           \
        }                                                               \
        else /* if (I64 && J64) */                                      \
        {                                                               \
            for (int64_t k = 0 ; k < nvals ; k++)                       \
            {                                                           \
                GrB_Index i = I64 [k] ;                                 \
                GrB_Index j = J64 [k] ;                                 \
                GRB_TRY (GrB_Matrix_setElement (A, XX [k], i, j)) ;     \
            }                                                           \
        }                                                               \
    }

    LG_SET_BURBLE (true) ;

    t = LAGraph_WallClockTime ( ) ;
    if (atype == GrB_BOOL)
    {
        BUILD_WITH_SETEL (bool) ;
    }
    else if (atype == GrB_INT8)
    {
        BUILD_WITH_SETEL (int8_t) ;
    }
    else if (atype == GrB_INT16)
    {
        BUILD_WITH_SETEL (int16_t) ;
    }
    else if (atype == GrB_INT32)
    {
        BUILD_WITH_SETEL (int32_t) ;
    }
    else if (atype == GrB_INT64)
    {
        BUILD_WITH_SETEL (int64_t) ;
    }
    else if (atype == GrB_UINT8)
    {
        BUILD_WITH_SETEL (uint8_t) ;
    }
    else if (atype == GrB_UINT16)
    {
        BUILD_WITH_SETEL (uint16_t) ;
    }
    else if (atype == GrB_UINT32)
    {
        BUILD_WITH_SETEL (uint32_t) ;
    }
    else if (atype == GrB_UINT64)
    {
        BUILD_WITH_SETEL (uint64_t) ;
    }
    else if (atype == GrB_FP32)
    {
        BUILD_WITH_SETEL (float) ;
    }
    else if (atype == GrB_FP64)
    {
        BUILD_WITH_SETEL (double) ;
    }
    else
    {
        printf ("type not supported\n") ;
        LG_ASSERT (false, GrB_INVALID_VALUE) ;
    }
    t = LAGraph_WallClockTime ( ) - t ;
    LG_SET_BURBLE (false) ;

    printf ("\n") ;
    printf ("Time for setElement:    %g sec\n", t) ;

    LG_SET_BURBLE (true) ;

    t = LAGraph_WallClockTime ( ) ;
    GRB_TRY (GrB_wait (A, GrB_MATERIALIZE)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for wait:          %g sec\n", t) ;

    LG_SET_BURBLE (false) ;

    //--------------------------------------------------------------------------
    // check the results (make sure A is a copy of G->A)
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    LG_TRY (LAGraph_Matrix_IsEqual (&isequal, A, G->A, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check:          %g sec\n", t) ;
    LG_ASSERT (isequal, GrB_INVALID_VALUE) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
#endif
    return (GrB_SUCCESS) ;
}

