//------------------------------------------------------------------------------
// mmtest: create a matrix and test Matrix Market I/O
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Michel Pelletier

#include "LAGraph.h"
#include <float.h>
#include <complex.h>
#include <assert.h>

#define OK(method)                                                      \
{                                                                       \
    GrB_Info this_info = method ;                                       \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))      \
    {                                                                   \
        printf ("complextest failure: [%d] %s\n", this_info, GrB_error ( )) ; \
        FREE_ALL ;                                                      \
        return (this_info) ;                                            \
    }                                                                   \
}                                                                       \

#define FREE_ALL                    \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&B) ;                 \
    GrB_free (&C) ;                 \
    GrB_free (&APPROX_ComplexFP64) ;    \
}                                   \

#define FOPENR(filename)                                        \
{                                                               \
    f = fopen (filename, "r") ;                                 \
    if (f == NULL)                                              \
    {                                                           \
        printf ("unable to open file [%s]\n", filename) ;       \
        FREE_ALL ;                                              \
        return (GrB_INVALID_VALUE) ;                            \
    }                                                           \
}                                                               \

#define CHECK(filename)                                     \
{                                                           \
    FOPENR(filename) ;                                      \
    OK (LAGraph_mmread (&D, f)) ;                           \
    fclose(f);                                              \
    OK (LAGraph_isequal(&bval, C, D, APPROX_ComplexFP64));      \
    if (!bval)                                              \
      exit(1);                                              \
}                                                           \

#define GET(M, I, J)                            \
  {                                             \
    OK(GrB_Matrix_extractElement_UDT(           \
    &val,                                       \
    M,                                          \
    I,                                          \
    J));                                        \
  }                                             \

#define GET_BOOL(M, I, J)                       \
  {                                             \
    OK(GrB_Matrix_extractElement_BOOL(          \
    &bval,                                      \
    M,                                          \
    I,                                          \
    J));                                        \
  }                                             \

#define GET_DOUBLE(M, I, J)                     \
  {                                             \
    OK(GrB_Matrix_extractElement_FP64(          \
    &dval,                                      \
    M,                                          \
    I,                                          \
    J));                                        \
  }                                             \

#define SET(M, I, J, V)                         \
  {                                             \
    val = V;                                    \
    OK(GrB_Matrix_setElement_UDT(               \
    M,                                          \
    &val,                                       \
    I,                                          \
    J));                                        \
  }                                             \

#define SET_DOUBLE(M, I, J, V)                  \
  {                                             \
    dval = V;                                   \
    OK(GrB_Matrix_setElement_FP64(              \
    M,                                          \
    dval,                                       \
    I,                                          \
    J));                                        \
  }                                             \

#define ZERO (CMPLX(0, 0))
#define ONE (CMPLX(1, 0))
#define LL (CMPLX(1.0, 1.0))
#define RR (CMPLX(2.0, 2.0))

#define TEST_BINOP(L, R, OP, V)                             \
  {                                                         \
    SET(A, 1, 1, L);                                        \
    SET(B, 1, 1, R);                                        \
    OK(GrB_eWiseAdd                                         \
       (C, NULL, NULL, OP, A, B, NULL)) ;                   \
    GET(C, 1, 1);                                           \
    printf("%s, %f + i%f\n", #OP, creal(val), cimag(val));  \
    if (val != V)                                           \
      exit(1);                                              \
  }                                                         \

#define TEST_BINOP_BOOL(L, R, OP, V)                        \
  {                                                         \
    SET(A, 1, 1, L);                                        \
    SET(B, 1, 1, R);                                        \
    OK(GrB_eWiseMult                                        \
       (C, NULL, NULL, OP, A, B, NULL)) ;                   \
    GET_BOOL(C, 1, 1);                                      \
    printf("%s, %i\n", #OP, bval);                          \
    if (bval != V)                                          \
      exit(1);                                              \
  }                                                         \

#define TEST_UOP(L, OP, V)                                  \
  {                                                         \
    SET(A, 1, 1, L);                                        \
    OK(GrB_apply                                            \
       (C, NULL, NULL, OP, A, NULL)) ;                      \
    GET(C, 1, 1);                                           \
    printf("%s, %f + i%f\n", #OP, creal(val), cimag(val));  \
    if (val != V)                                           \
      exit(1);                                              \
  }                                                         \

#define TEST_UOP_BOOL(L, OP, V)                             \
  {                                                         \
    SET(A, 1, 1, L);                                        \
    OK(GrB_apply                                            \
       (C, NULL, NULL, OP, A, NULL)) ;                      \
    GET_BOOL(C, 1, 1);                                      \
    printf("%s, %i\n", #OP, bval);                          \
    if (bval != V)                                          \
      exit(1);                                              \
  }                                                         \

#define TEST_UOP_DOUBLE(L, OP, V)                           \
  {                                                         \
    SET(A, 1, 1, L);                                        \
    OK(GrB_apply                                            \
       (C, NULL, NULL, OP, A, NULL)) ;                      \
    GET_DOUBLE(C, 1, 1);                                    \
    printf("%s, %f\n", #OP, dval);                          \
    if (dval != V)                                          \
      exit(1);                                              \
  }                                                         \

#define TEST_UOP_CMPLX(L, OP, V)                            \
  {                                                         \
    SET_DOUBLE(A, 1, 1, L);                                 \
    OK(GrB_apply                                            \
       (C, NULL, NULL, OP, A, NULL)) ;                      \
    GET(C, 1, 1);                                           \
    printf("%s, %f + i%f\n", #OP, creal(val), cimag(val));  \
    if (val != V)                                           \
      exit(1);                                              \
  }                                                         \

void complex_approx
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    *z = cabs(*x - *y) <= DBL_EPSILON;
}

//------------------------------------------------------------------------------
// complextest main:
//------------------------------------------------------------------------------

int main (int argc, char **argv)
{

    printf ("Complex/complex: "
        "testing LAGraph_ComplexFP64 and its operators:\n") ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    printf ("LAGraph_xinit (requires SuiteSparse:GraphBLAS)\n") ;
    LAGraph_xinit (malloc, calloc, realloc, free, true) ;
    #else
    printf ("LAGraph_init\n") ;
    LAGraph_init ( ) ;
    #endif

    GrB_BinaryOp APPROX_ComplexFP64;

    GrB_Matrix A = NULL, B = NULL, C = NULL, D = NULL, E = NULL;
    uint64_t aseed = 42;
    uint64_t bseed = 43;
    double complex val = 3.0 + 4.0i;
    bool bval = false;
    double dval = 0.0;
    FILE *f = NULL ;

    OK (GrB_BinaryOp_new(
         &APPROX_ComplexFP64,
         (LAGraph_binary_function) (&complex_approx),
         GrB_BOOL, LAGraph_ComplexFP64, LAGraph_ComplexFP64)) ;

    OK(LAGraph_random
       (&A,
        LAGraph_ComplexFP64,
        2,
        2,
        3,
        false,
        false,
        false,
        false,
        false,
        &aseed
        )) ;

    OK(LAGraph_random
       (&B,
        LAGraph_ComplexFP64,
        2,
        2,
        3,
        false,
        false,
        false,
        false,
        false,
        &bseed
        )) ;

    OK(GrB_Matrix_new
       (&C, LAGraph_ComplexFP64, 2, 2)) ;

    OK(GrB_Matrix_new
       (&E, GrB_BOOL, 2, 2)) ;

    // Add matrices

    OK(GrB_eWiseAdd
       (C, NULL, NULL, LAGraph_PLUS_ComplexFP64, A, B, NULL)) ;

    CHECK("data/test_eadd.mtx");

    /* // Multi matrices */
    OK(GrB_eWiseMult
       (C, NULL, NULL, LAGraph_TIMES_ComplexFP64, A, B, NULL)) ;

    CHECK("data/test_emul.mtx");

    // Matrix Mult matrices
    OK(GrB_mxm
       (C, NULL, NULL, LAGraph_PLUS_TIMES_ComplexFP64, A, B, NULL)) ;

    CHECK("data/test_mxm.mtx");

    TEST_BINOP(LL, RR, LAGraph_MAX_ComplexFP64, RR);
    TEST_BINOP(LL, RR, LAGraph_MIN_ComplexFP64, LL);
    TEST_BINOP(LL, RR, LAGraph_FIRST_ComplexFP64, LL);
    TEST_BINOP(LL, RR, LAGraph_SECOND_ComplexFP64, RR);
    TEST_BINOP(LL, RR, LAGraph_PLUS_ComplexFP64, CMPLX(3, 3));
    TEST_BINOP(LL, RR, LAGraph_MINUS_ComplexFP64, CMPLX(-1, -1));
    TEST_BINOP(LL, RR, LAGraph_RMINUS_ComplexFP64, LL);
    TEST_BINOP(LL, RR, LAGraph_TIMES_ComplexFP64, CMPLX(0, 4));
    TEST_BINOP(LL, RR, LAGraph_DIV_ComplexFP64, CMPLX(0.5, 0));
    TEST_BINOP(LL, RR, LAGraph_RDIV_ComplexFP64, CMPLX(2, 0));
    TEST_BINOP(LL, RR, LAGraph_PAIR_ComplexFP64, CMPLX(1, 0));
    TEST_BINOP(LL, RR, LAGraph_ANY_ComplexFP64, RR);
    TEST_BINOP(LL, RR, LAGraph_ISEQ_ComplexFP64, ZERO);
    TEST_BINOP(LL, RR, LAGraph_ISNE_ComplexFP64, ONE);
    TEST_BINOP(LL, RR, LAGraph_ISGT_ComplexFP64, ZERO);
    TEST_BINOP(LL, RR, LAGraph_ISLT_ComplexFP64, ONE);
    TEST_BINOP(LL, RR, LAGraph_ISGE_ComplexFP64, ZERO);
    TEST_BINOP(LL, RR, LAGraph_ISLE_ComplexFP64, ONE);
    TEST_BINOP(LL, RR, LAGraph_OR_ComplexFP64, ONE);
    TEST_BINOP(LL, RR, LAGraph_AND_ComplexFP64, ONE);
    TEST_BINOP(LL, RR, LAGraph_XOR_ComplexFP64, ZERO);

    TEST_UOP(LL, LAGraph_ONE_ComplexFP64, ONE);
    TEST_UOP(RR, LAGraph_IDENTITY_ComplexFP64, RR);
    TEST_UOP(RR, LAGraph_AINV_ComplexFP64, CMPLX(-2, -2));
    TEST_UOP(CMPLX(-2, 0), LAGraph_ABS_ComplexFP64, CMPLX(2, 0));
    TEST_UOP(CMPLX(-2, 0), LAGraph_MINV_ComplexFP64, CMPLX(-0.5, -0));
    TEST_UOP(CMPLX(-2, 0), LAGraph_NOT_ComplexFP64, CMPLX(0, 0));
    TEST_UOP(CMPLX(-2, 2), LAGraph_CONJ_ComplexFP64, CMPLX(-2, -2));

    OK(GrB_free(&C));
    OK(GrB_Matrix_new
       (&C, GrB_BOOL, 2, 2)) ;

    TEST_BINOP_BOOL(LL, RR, LAGraph_EQ_ComplexFP64, false);
    TEST_BINOP_BOOL(LL, RR, LAGraph_NE_ComplexFP64, true);
    TEST_BINOP_BOOL(LL, RR, LAGraph_GT_ComplexFP64, false);
    TEST_BINOP_BOOL(LL, RR, LAGraph_LT_ComplexFP64, true);
    TEST_BINOP_BOOL(LL, RR, LAGraph_GE_ComplexFP64, false);
    TEST_BINOP_BOOL(LL, RR, LAGraph_LE_ComplexFP64, true);
    TEST_BINOP_BOOL(LL, RR, LAGraph_SKEW_ComplexFP64, false);
    TEST_BINOP_BOOL(LL, RR, LAGraph_HERMITIAN_ComplexFP64, false);
    TEST_UOP_BOOL(CMPLX(1, 0), LAGraph_ISONE_ComplexFP64, true);
    TEST_UOP_BOOL(CMPLX(-2, 2), LAGraph_TRUE_BOOL_ComplexFP64, true);

    OK(GrB_free(&C));
    OK(GrB_Matrix_new
       (&C, GrB_FP64, 2, 2)) ;

    TEST_UOP_DOUBLE(CMPLX(-2, 0), LAGraph_REAL_ComplexFP64, -2);
    TEST_UOP_DOUBLE(CMPLX(-2, 2), LAGraph_IMAG_ComplexFP64, 2);
    TEST_UOP_DOUBLE(CMPLX(-2, 0), LAGraph_CABS_ComplexFP64, 2);
    TEST_UOP_DOUBLE(CMPLX(1, 0), LAGraph_ANGLE_ComplexFP64, 0);

    OK(GrB_free(&A));
    OK(GrB_Matrix_new
       (&A, GrB_FP64, 2, 2)) ;

    OK(GrB_free(&C));
    OK(GrB_Matrix_new
       (&C, LAGraph_ComplexFP64, 2, 2)) ;

    TEST_UOP_CMPLX(-2, LAGraph_COMPLEX_REAL_ComplexFP64, CMPLX(-2, 0));
    TEST_UOP_CMPLX(2, LAGraph_COMPLEX_IMAG_ComplexFP64, CMPLX(0, 2));
}
