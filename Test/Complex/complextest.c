#include "LAGraph.h"
#include <float.h>
#include <complex.h>
#include <assert.h>

#define OK(method) do                                                   \
{                                                                       \
    GrB_Info this_info = method ;                                       \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))      \
    {                                                                   \
        printf ("complextest failure: [%d] %s\n", this_info, GrB_error ( )) ; \
        FREE_ALL ;                                                      \
        return (this_info) ;                                            \
    }                                                                   \
} while(0);

#define FREE_ALL do                 \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&B) ;                 \
    GrB_free (&C) ;                 \
    GrB_free (&APPROX_Complex) ;    \
} while(0);

#define FOPENR(filename) do                                     \
{                                                               \
    f = fopen (filename, "r") ;                                 \
    if (f == NULL)                                              \
    {                                                           \
        printf ("unable to open file [%s]\n", filename) ;       \
        FREE_ALL ;                                              \
        return (GrB_INVALID_VALUE) ;                            \
    }                                                           \
} while(0) ;

#define CHECK(filename) do {                        \
    FOPENR(filename) ;                                      \
    OK (LAGraph_mmread (&D, f)) ;                          \
    fclose(f);                                              \
    OK (LAGraph_isequal(&result, C, D, APPROX_Complex));  \
    assert(result);                                         \
} while(0) ;


//------------------------------------------------------------------------------
// mmtest main:
//------------------------------------------------------------------------------

// calls mmtest with a huge range of combinations of matrix sizes, types,
// and characteristics.

void complex_approx
(
    bool *z,
    const double complex *x,
    const double complex *y
)
{
    *z = cabs(*x - *y) <= DBL_EPSILON;
}


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
