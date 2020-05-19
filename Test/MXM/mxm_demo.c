//------------------------------------------------------------------------------
// GraphBLAS/alternative/mxm_demo.c: test GrB_mxm
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2020, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Usage:
//
//  mxm_demo

#include "LAGraph.h"
#include "/home/davis/sparse/GraphBLAS/Source/GB_Global.h"

#define LAGRAPH_FREE_ALL    \
    LAGr_free (&Cout) ;     \
    LAGr_free (&Cin) ;      \
    LAGr_free (&M) ;        \
    LAGr_free (&A) ;        \
    LAGr_free (&B) ;

int main (int argc, char **argv)
{
    GrB_Matrix Cin = NULL, Cout = NULL, A = NULL, B = NULL, M = NULL ;
    GrB_Info info ;
    double tic [2], r1, r2 ;
    LAGraph_init () ;
    int nthreads ;
    LAGr_get (GxB_NTHREADS, &nthreads) ;
    fprintf (stderr, "mxm_demo: nthreads %d\n", nthreads) ;
    printf ("--------------------------------------------------------------\n");

    // GraphBLAS semirings supported by MKL
    #define NSEMIRINGS 13
    GrB_Semiring Semirings [NSEMIRINGS] =
    {
        GrB_LOR_LAND_SEMIRING_BOOL,
        GrB_PLUS_TIMES_SEMIRING_INT32,
        GrB_PLUS_TIMES_SEMIRING_INT64,
        GrB_PLUS_TIMES_SEMIRING_FP32,
        GrB_PLUS_TIMES_SEMIRING_FP64,
        GrB_MIN_PLUS_SEMIRING_INT32,
        GrB_MIN_PLUS_SEMIRING_INT64,
        GrB_MIN_PLUS_SEMIRING_FP32,
        GrB_MIN_PLUS_SEMIRING_FP64,
        GrB_MAX_FIRST_SEMIRING_INT32,
        GrB_MAX_FIRST_SEMIRING_INT64,
        GrB_MAX_FIRST_SEMIRING_FP32,
        GrB_MAX_FIRST_SEMIRING_FP64,
    } ;

    // corresponding data types for the 13 semirings:
    GrB_Type Types [NSEMIRINGS] =
    {
        GrB_BOOL,
        GrB_INT32, GrB_INT64, GrB_FP32, GrB_FP64,
        GrB_INT32, GrB_INT64, GrB_FP32, GrB_FP64,
        GrB_INT32, GrB_INT64, GrB_FP32, GrB_FP64,
    } ;

    // C = A*B where A is m-by-k and B is k-by-n
    #define NPROBLEMS 5
    int64_t m [NPROBLEMS] = { 1, 2, 10, 4, 5  } ;
    int64_t n [NPROBLEMS] = { 1, 3, 10, 2, 4  } ;
    int64_t k [NPROBLEMS] = { 1, 5, 10, 1, 6  } ;
    int64_t v [NPROBLEMS] = { 3, 7, 30, 5, 20 } ;

    #define RANDOM(A,type,m,n,v) \
    {                                                                       \
        LAGr_random (&A, type, m, n, v, false, false, false, false, false,  \
            &seed) ;                                                        \
        GrB_Index ignore  ;                                                 \
        GrB_Matrix_nvals (&ignore, A) ;                                     \
    }

    //--------------------------------------------------------------------------
    // test each semiring
    //--------------------------------------------------------------------------

    for (int s = 0 ; s < NSEMIRINGS ; s++)
    {
        GrB_Semiring semiring = Semirings [s] ;
        GrB_Type type = Types [s] ;

        printf ("\n======================================================\n") ;
        GxB_print (semiring, 3) ;
        GxB_print (type, 3) ;
        printf ("\n======================================================\n") ;

        //----------------------------------------------------------------------
        // test each problem
        //----------------------------------------------------------------------

        for (int t = 0 ; t < NPROBLEMS ; t++)
        {

            //------------------------------------------------------------------
            // create Cin, M, A, and B
            //------------------------------------------------------------------

            printf ("\n    ----------------------------------------------\n") ;
            uint64_t seed = t ;
            RANDOM (Cin, type, m [t], n [t], v [t]) ;
            RANDOM (M,   type, m [t], n [t], v [t]) ;
            RANDOM (A,   type, m [t], k [t], v [t]) ;
            RANDOM (B,   type, k [t], n [t], v [t]) ;

            GxB_print (Cin, 3) ;
            GxB_print (M, 3) ;
            GxB_print (A, 3) ;
            GxB_print (B, 3) ;

            //------------------------------------------------------------------
            // C = A*B
            //------------------------------------------------------------------

            LAGr_Matrix_dup (&Cout, Cin) ;
            LAGr_mxm (Cout, NULL, NULL, semiring, A, B, NULL) ;
            GxB_print (Cout, 3) ;
            LAGr_free (&Cout) ;

            LAGRAPH_FREE_ALL  ;
        }
    }
}

