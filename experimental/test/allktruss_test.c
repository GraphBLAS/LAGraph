//------------------------------------------------------------------------------
// LAGraph/experimental/test/allktruss_test.c: test program for LAGraph_ktruss
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Usage:  allktest < matrixmarketfile.mtx
// Contributed by Tim Davis, Texas A&M

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL                            \
{                                                   \
    if (Cset != NULL)                               \
    {                                               \
        for (int64_t kk = 3 ; kk <= kmax ; kk++)    \
        {                                           \
            GrB_free (&(Cset [kk])) ;               \
        }                                           \
    }                                               \
    GrB_free (&A) ;                                 \
    GrB_free (&M) ;                                 \
    GrB_free (&LAGraph_ONE_UINT32) ;                \
    GrB_free (&LAGraph_LOR_UINT32) ;                \
    LAGraph_Free ((void **) &ntris) ;                          \
    LAGraph_Free ((void **) &nedges) ;                         \
    LAGraph_Free ((void **) &nstepss) ;                        \
}

//****************************************************************************

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

void LAGraph_one_uint32
(
    uint32_t *z,
    const uint32_t *x       // ignored
)
{
    (*z) = 1 ;
}

void LAGraph_lor_uint32
(
    uint32_t *z,
    const uint32_t *x,
    const uint32_t *y
)
{
    (*z) = (((*x) != 0) || ((*y) != 0)) ;
}

//****************************************************************************
int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, C = NULL, M = NULL ;
    GrB_Type C_type = NULL;
    GrB_Matrix *Cset = NULL ;
    GrB_UnaryOp  LAGraph_ONE_UINT32 = NULL;
    GrB_BinaryOp LAGraph_LOR_UINT32 = NULL;
    int64_t kmax = 0 ;
    int64_t *ntris = NULL, *nedges  = NULL, *nstepss = NULL ;

    LAGRAPH_OK (LAGraph_Init (NULL)) ;
    int nthreads_max;
    LAGRAPH_OK (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;

#if 0 //LG_SUITESPARSE
    LAGraph_ONE_UINT32 = GxB_ONE_UINT32 ;
    LAGraph_LOR_UINT32 = GxB_LOR_UINT32 ;
#else
    // create a new built-in unary operator using LAGraph_one_uint32
    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ONE_UINT32,
        F_UNARY (LAGraph_one_uint32),
        GrB_UINT32, GrB_UINT32)) ;

    // create a new built-in binary operator using LAGraph_lor_uint32
    LAGRAPH_OK (GrB_BinaryOp_new (&LAGraph_LOR_UINT32,
        F_BINARY (LAGraph_lor_uint32),
        GrB_UINT32, GrB_UINT32, GrB_UINT32)) ;
#endif

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGRAPH_OK (LAGraph_Tic (tic, NULL)) ;

    FILE *f ;
    if (argc == 1)
    {
        f = stdin ;
    }
    else
    {
        f = fopen (argv[1], "r") ;
        if (f == NULL)
        {
            printf ("unable to open file [%s]\n", argv[1]) ;
            return (GrB_INVALID_VALUE) ;
        }
    }
    LAGRAPH_OK (LAGraph_MMRead (&C, &C_type, f, NULL)) ;
    double t_read;
    LAGRAPH_OK (LAGraph_Toc (&t_read, tic, NULL)) ;
    printf ("\nread A time:     %14.6f sec\n", t_read) ;

    LAGRAPH_OK (LAGraph_Tic (tic, NULL)) ;
    GrB_Index n, ne ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, C)) ;

    // A = spones (C), and typecast to uint32
    LAGRAPH_OK (GrB_Matrix_new (&A, GrB_UINT32, n, n)) ;
    LAGRAPH_OK (GrB_apply (A, NULL, NULL, LAGraph_ONE_UINT32, C, NULL)) ;
    GrB_free (&C) ;

    // M = diagonal mask matrix
    LAGRAPH_OK (GrB_Matrix_new (&M, GrB_BOOL, n, n)) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        // M(i,i) = true ;
        LAGRAPH_OK (GrB_Matrix_setElement (M, (bool) true, i, i)) ;
    }

    // make A symmetric (A = spones (A+A')) and remove self edges (via M)
    LAGRAPH_OK (GrB_eWiseAdd (A, M, NULL, LAGraph_LOR_UINT32, A, A,
        GrB_DESC_RCT1)) ;
    GrB_free (&M) ;

    LAGRAPH_OK (GrB_Matrix_nvals (&ne, A)) ;

    double t_process;
    LAGRAPH_OK (LAGraph_Toc (&t_process, tic, NULL)) ;
    printf ("process A time:  %14.6f sec\n", t_process) ;

    //--------------------------------------------------------------------------
    // construct all ktrusses
    //--------------------------------------------------------------------------

    Cset = NULL ; // LAGraph_malloc (n, sizeof (GrB_Matrix)) ;
    ntris   = LAGraph_Malloc (n, sizeof (int64_t)) ;
    nedges  = LAGraph_Malloc (n, sizeof (int64_t)) ;
    nstepss = LAGraph_Malloc (n, sizeof (int64_t)) ;

    if (ntris == NULL || nedges == NULL || nstepss == NULL)
    {
        printf ("out of memory\n") ;
        LAGraph_FREE_ALL ;
    }

    double t1 ;
    for (int nthreads = 1 ; nthreads <= nthreads_max ; )
    {
        LAGRAPH_OK (LAGraph_SetNumThreads (nthreads, NULL)) ;

        double tic [2] ;
        LAGRAPH_OK (LAGraph_Tic (tic, NULL)) ;
        LAGRAPH_OK (LAGraph_allktruss (Cset, A, &kmax, ntris, nedges, nstepss));
        double t;
        LAGRAPH_OK (LAGraph_Toc (&t, tic, NULL)) ;

        if (nthreads == 1)
        {
            t1 = t ;
            for (int64_t kk = 3 ; kk <= kmax ; kk++)
            {
                printf (" k %4d", (int) kk) ;
                printf (" edges %12" PRIu64, nedges [kk]) ;
                printf (" ntriangles %12" PRId64, ntris [kk]) ;
                printf (" nsteps %6" PRId64 "\n", nstepss [kk]) ;
            }
        }

        printf ("nthreads: %3d time: %12.6f rate: %6.2f", nthreads, t,
            1e-6 * ne / t) ;
        if (nthreads > 1)
        {
            printf (" speedup: %6.2f", t1 / t) ;
        }
        printf ("\n") ;

        if (nthreads != nthreads_max && 4 * nthreads > nthreads_max)
        {
            nthreads = nthreads_max ;
        }
        else
        {
            nthreads *= 4 ;
        }
    }

    printf ("\n") ;
    LAGraph_FREE_ALL ;
    LAGraph_Finalize (NULL) ;
}
