#include "LAGraph.h"

#define OK(method)                                                          \
{                                                                           \
  printf(#method);                                                     \
  printf("\n");                                                         \
    GrB_Info this_info = method ;                                           \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))          \
    {                                                                       \
        printf ("complextest failure: [%d] %s\n", this_info, GrB_error ( )) ;    \
        FREE_ALL ;                                                          \
        return (this_info) ;                                                \
    }                                                                       \
}

#define FREE_ALL                    \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&B) ;                 \
}

//------------------------------------------------------------------------------
// mmtest main:
//------------------------------------------------------------------------------

// calls mmtest with a huge range of combinations of matrix sizes, types,
// and characteristics.

int main (int argc, char **argv)
{

    printf ("Complex/complex: "
        "testing LAGraph_Complex and its operators:\n") ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    printf ("testing LAGraph_xinit (requires SuiteSparse:GraphBLAS)\n") ;
    LAGraph_xinit (malloc, calloc, realloc, free, true) ;
    #else
    printf ("LAGraph_init\n") ;
    LAGraph_init ( ) ;
    #endif

    // Create random complex matrix

    GrB_Matrix A = NULL, B = NULL, C = NULL;
    uint64_t aseed = 42;
    uint64_t bseed = 43;
    double complex val = 3.0 + 4.0i;

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
        ));

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
        ));

    OK(GrB_Matrix_new
       (&C, LAGraph_Complex, 2, 2));

    OK (LAGraph_mmwrite (A, stdout)) ;
    OK (LAGraph_mmwrite (B, stdout)) ;
    
    // Add matrices

    OK(GrB_eWiseAdd
       (C, NULL, NULL, LAGraph_PLUS_Complex, A, B, NULL));

    OK (LAGraph_mmwrite (C, stdout)) ;
    
    /* // Multi matrices */
    /* OK(GrB_eWiseMult */
    /*    (C, NULL, NULL, LAGraph_TIMES_Complex, A, B, NULL)); */

    /* OK (LAGraph_mmwrite (C, stdout)) ; */
    
    // Matrix Mult matrices
    OK(GrB_mxm
       (C, NULL, NULL, LAGraph_PLUS_TIMES_Complex, A, B, NULL));

    OK (LAGraph_mmwrite (C, stdout)) ;
    
}
