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
    GrB_free (&APPROX_Complex) ;    \
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
    OK (LAGraph_isequal(&result, C, D, APPROX_Complex));    \
    assert(result);                                         \
}                                                           \


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
        "testing LAGraph_Complex and its operators:\n") ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    printf ("LAGraph_xinit (requires SuiteSparse:GraphBLAS)\n") ;
    LAGraph_xinit (malloc, calloc, realloc, free, true) ;
    #else
    printf ("LAGraph_init\n") ;
    LAGraph_init ( ) ;
    #endif

    GrB_BinaryOp APPROX_Complex;

    GrB_Matrix A = NULL, B = NULL, C = NULL, D = NULL;
    uint64_t aseed = 42;
    uint64_t bseed = 43;
    double complex val = 3.0 + 4.0i;
    bool result = false;
    FILE *f = NULL ;

    OK (GrB_BinaryOp_new(
         &APPROX_Complex,
         (LAGraph_binary_function) (&complex_approx),
         GrB_BOOL, LAGraph_Complex, LAGraph_Complex)) ;

    OK(LAGraph_random
       (&A,
        LAGraph_Complex,
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
        LAGraph_Complex,
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
       (&C, LAGraph_Complex, 2, 2)) ;

    // Add matrices

    OK(GrB_eWiseAdd
       (C, NULL, NULL, LAGraph_PLUS_Complex, A, B, NULL)) ;

    CHECK("data/test_eadd.mtx");

    /* // Multi matrices */
    OK(GrB_eWiseMult
       (C, NULL, NULL, LAGraph_TIMES_Complex, A, B, NULL)) ;

    CHECK("data/test_emul.mtx");

    // Matrix Mult matrices
    OK(GrB_mxm
       (C, NULL, NULL, LAGraph_PLUS_TIMES_Complex, A, B, NULL)) ;

    CHECK("data/test_mxm.mtx");
}
