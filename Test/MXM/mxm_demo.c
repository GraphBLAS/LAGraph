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
#include "../../../GraphBLAS/Source/GB_Global.h"

#define LAGRAPH_FREE_ALL    \
    (*ok) = false ;         \
    FREE_WORK ;

#define FREE_WORK           \
    LAGr_free (&S1) ;       \
    LAGr_free (&S2) ;       \
    LAGr_free (&E) ;

GrB_Info check_results
(
    bool *ok,
    GrB_Matrix C1,
    GrB_Matrix C2,
    GrB_Type type
)
{
    GrB_Info info ;
    GrB_Matrix S1 = NULL, S2 = NULL, E = NULL ;
    (*ok) = false ;
    if (type == GrB_FP32 || type == GrB_FP64)
    {
        // pattern must be equal.  values can differ by O(eps)
        LAGRAPH_OK (LAGraph_pattern (&S1, C1, NULL)) ;
        LAGRAPH_OK (LAGraph_pattern (&S2, C2, NULL)) ;
        LAGRAPH_OK (LAGraph_isequal (ok, S1, S2, NULL)) ;
        LAGr_free (&S1) ;
        LAGr_free (&S2) ;
        if (!(*ok))
        {
            printf ("ERROR: pattern of C1 and C2 differ!\n") ;
        }
        // err = sum (abs (C1-C2))
        (*ok) = false ;
        GrB_Index m, n ;
        LAGr_Matrix_nrows (&m, C1) ;
        LAGr_Matrix_ncols (&n, C1) ;
        LAGr_Matrix_new (&E, GrB_FP64, m, n) ;
        LAGr_eWiseMult (E, NULL, NULL, GrB_MINUS_FP64, C1, C2, NULL) ;
        LAGr_apply (E, NULL, NULL, GrB_ABS_FP64, E, NULL) ;
        double err = 1 ;
        LAGr_reduce (&err, NULL, GrB_PLUS_MONOID_FP64, E, NULL) ;
        LAGr_free (&E) ;
        printf ("norm (C1-C2) = %g  ", err) ;
        (*ok) = (err < ((type == GrB_FP32) ? 1e-6 : 1e-12)) ;
        if (!(*ok))
        {
            printf ("ERROR: norm too high!\n") ;
        }
        printf ("\n") ;
    }
    else
    {
        // matrices must be identical
        LAGRAPH_OK (LAGraph_isequal (ok, C1, C2, NULL)) ;
        if (!(*ok))
        {
            printf ("ERROR: C1 and C2 differ!\n") ;
        }
    }
    return (GrB_SUCCESS) ;
}

#undef  FREE_WORK
#define FREE_WORK           \
    LAGr_free (&C1) ;       \
    LAGr_free (&C2) ;       \
    LAGr_free (&Cin) ;      \
    LAGr_free (&M) ;        \
    LAGr_free (&A) ;        \
    LAGr_free (&B) ;

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL    \
    FREE_WORK ;             \
    LAGraph_finalize ( ) ;

int main (int argc, char **argv)
{
    GrB_Matrix Cin = NULL, C1 = NULL, A = NULL, B = NULL, M = NULL, C2 = NULL ;
    GrB_Info info ;
    double tic [2], r1, r2 ;
    LAGraph_init () ;
    GxB_set (GxB_BURBLE, true) ;
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

    int nfail = 0 ;

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

            for (int mask = 0 ; mask <= 1 ; mask++)
            {
                printf ("mask:%d\n", mask) ;
                GrB_Matrix M1 = mask ? M : NULL ;

                //--------------------------------------------------------------
                // C1=A*B or C1<M>=A*B without MKL
                //--------------------------------------------------------------

                GB_Global_hack_set (0) ;

                LAGr_Matrix_dup (&C1, Cin) ;
                LAGr_mxm (C1, M1, NULL, semiring, A, B, NULL) ;
                GxB_print (C1, 3) ;

                //--------------------------------------------------------------
                // C2=A*B or C2<M>=A*B without MKL
                //--------------------------------------------------------------

                GB_Global_hack_set (1) ;

                LAGr_Matrix_dup (&C2, Cin) ;
                LAGr_mxm (C2, M1, NULL, semiring, A, B, NULL) ;
                GxB_print (C2, 3) ;

                //--------------------------------------------------------------
                // compare results
                //--------------------------------------------------------------

                bool ok = false ;
                LAGRAPH_OK (check_results (&ok, C1, C2, type))
                if (!ok) nfail++ ;
                LAGr_free (&C1) ;
                LAGr_free (&C2) ;
            }

            //------------------------------------------------------------------
            // free all matrices
            //------------------------------------------------------------------

            FREE_WORK  ;
        }
    }

    printf ("test failures: %d\n", nfail) ;
    LAGraph_finalize ( ) ;
}

