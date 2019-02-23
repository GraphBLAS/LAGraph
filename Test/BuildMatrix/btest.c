//------------------------------------------------------------------------------
// btest: create a matrix with lots of threads
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

// This test requires SuiteSparse:GraphBLAS

//------------------------------------------------------------------------------

#include "LAGraph.h"

#define OK(method)                                                          \
{                                                                           \
    GrB_Info this_info = method ;                                           \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))          \
    {                                                                       \
        printf ("btest failure: [%d] %s\n", this_info, GrB_error ( )) ;     \
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
// bmake: create a matrix with lots of threads
//------------------------------------------------------------------------------

GrB_Info bmake
(
    GrB_Matrix *A_handle,       // matrix created
    GrB_Type type,
    GrB_Index nrows,
    GrB_Index ncols,
    GrB_Index nvals,
    bool make_pattern,
    bool make_symmetric,
    bool make_skew_symmetric,
    bool make_hermitian,
    bool no_diagonal,
    uint64_t *seed,
    int nthreads                // number of threads to use
)
{
    GrB_Matrix A = NULL, B = NULL ;

    //--------------------------------------------------------------------------
    // create a random matrix
    //--------------------------------------------------------------------------

    OK (GxB_set (GxB_NTHREADS, nthreads)) ;

    // this is with one thread, so far:
    OK (LAGraph_random (&A, type, nrows, ncols, nvals,
        make_pattern, make_symmetric, make_skew_symmetric,
        make_hermitian, no_diagonal, seed)) ;

    // finish any pending computations (with nthreads in the qsort)
    // printf ("A random before sort:\n") ;
    // OK (GxB_print (A, 2)) ;

    double tic [2], t ;
    LAGraph_tic (tic) ;

    OK (GrB_Matrix_nvals (&nvals, A)) ;

    t = LAGraph_toc (tic) ;
    printf ("A random after sort:  nthreads %d time %g\n", nthreads, t) ;
    // OK (GxB_print (A, 2)) ;

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*A_handle) = A ;
    return (GrB_SUCCESS) ;
}


//------------------------------------------------------------------------------
// btest main:
//------------------------------------------------------------------------------

int main (int argc, char **argv)
{

    GrB_Matrix A = NULL, B = NULL ;

    printf ("BuildMatrix/btest: test LAGraph_random with parallel qsort\n") ;

    LAGraph_init ( ) ;

    // TODO: put this in the LAGraph_Context, created by LAGraph_init
    #ifdef _OPENMP
    int maxthreads = omp_get_max_threads ( ) ;
    #else
    int maxthreads = 1 ;
    #endif

    printf ("max threads %d\n", maxthreads) ;

    uint64_t seed = 1 ;
    GrB_Index nrows = 100000 ;
    GrB_Index ncols = 100000 ;
    GrB_Index nvals = 10000000 ;

    // create A with one thread
    OK (bmake (&A, GrB_FP64, nrows, ncols, nvals,
        false, false, false, false, false, &seed, 1)) ;

    // OK (GxB_print (A, 2)) ;

    // now create B with 1 to max # of threads
    for (int nthreads = 1 ; nthreads <= maxthreads ; nthreads++)
    {

        // reset the seed
        seed = 1 ;

        // create B with nthreads
        OK (bmake (&B, GrB_FP64, nrows, ncols, nvals,
            false, false, false, false, false, &seed, nthreads)) ;

        // OK (GxB_print (B, 2)) ;

        bool ok = false ;
        OK (LAGraph_isequal (&ok, A, B, NULL)) ;

        if (!ok)
        {
            printf ("A and B do not match\n") ;
            // printf ("--------------A:\n") ;
            // GxB_print (A, GxB_COMPLETE) ;
            // OK (LAGraph_mmwrite (A, stdout)) ;
            // printf ("--------------B:\n") ;
            // GxB_print (B, GxB_COMPLETE) ;
            // OK (LAGraph_mmwrite (B, stdout)) ;
            FREE_ALL ;
            return (GrB_INVALID_VALUE) ;
        }

        OK (GrB_free (&B)) ;
    }

    printf ("\nbtest: all tests passed\n") ;
    FREE_ALL ;
    LAGraph_finalize ( ) ;
    return (GrB_SUCCESS) ;
}

