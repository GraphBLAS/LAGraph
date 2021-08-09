//------------------------------------------------------------------------------
// LAGraph/Test/KTruss/ktest.c: test program for LAGraph_ktruss
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage:  ktruss_test < matrixmarketfile.mtx
// Usage:  ktruss_test matrixmarketfile.mtx

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&C) ;         \
    GrB_free (&A) ;         \
    GrB_free (&M) ;         \
    GrB_free (&LAGraph_ONE_UINT32) ;                \
    GrB_free (&LAGraph_LOR_UINT32) ;                \
}

//****************************************************************************

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

void LAGraph_one_uint32
(
    double *z,
    const double *x       // ignored
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
    GrB_UnaryOp  LAGraph_ONE_UINT32 = NULL;
    GrB_BinaryOp LAGraph_LOR_UINT32 = NULL;

    LAGraph_Init (NULL) ;
    int nthreads_max;
    LAGraph_GetNumThreads (&nthreads_max, NULL) ;

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_Tic (tic, NULL) ;

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
    LAGraph_Toc (&t_read, tic, NULL) ;
    printf ("\nread A time:     %14.6f sec\n", t_read) ;

    LAGraph_Tic (tic, NULL) ;
    GrB_Index n, nedges ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, C)) ;

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
    GrB_free (&M) ; M = NULL;

    LAGRAPH_OK (GrB_Matrix_nvals (&nedges, A)) ;

    double t_process;
    LAGraph_Toc (&t_process, tic, NULL) ;
    printf ("process A time:  %14.6f sec\n", t_process) ;
    printf ("input graph: %g nodes, %g edges\n", (double) n, (double) nedges) ;

    //--------------------------------------------------------------------------
    // construct all ktrusses
    //--------------------------------------------------------------------------

    int64_t nt ;
    GrB_Index nedges_in_ktruss = 1 ;

    for (uint32_t k = 3 ; k < 10 && nedges_in_ktruss > 0 ; k++)
    {

        printf ("\nKTruss: k = %3d:\n", (int) k) ;
        double t1 ;

        double t_sequential ;
        for (int nthreads = 1 ; nthreads <= nthreads_max ; )
        {
            LAGraph_SetNumThreads (nthreads, NULL) ;

            double tic [2] ;
            LAGraph_Tic (tic, NULL) ;

            int32_t nsteps ;
            LAGRAPH_OK (LAGraph_ktruss (&C, &C_type, A, k, &nsteps)) ;
            LAGRAPH_OK (GrB_Matrix_nvals (&nedges_in_ktruss, C)) ;
            double t;
            LAGraph_Toc (&t, tic, NULL) ;

            if (nthreads == 1)
            {
                t1 = t ;
                LAGRAPH_OK (GrB_reduce (&nt, NULL, GrB_PLUS_MONOID_INT64,
                                        C, NULL)) ;
                nt /= 6 ;
                printf (" edges %" PRIu64, nedges_in_ktruss/2) ;
                printf (" ntriangles %" PRId64 "\n", nt) ;
            }

            GrB_free (&C) ; C = NULL;

            printf ("nthreads: %3d time: %12.6f rate: %6.2f", nthreads, t,
                    1e-6 * nedges / t) ;
            if (nthreads > 1)
            {
                printf (" speedup: %6.2f", t1 / t) ;
            }
            printf (" steps %d", nsteps) ;
            printf ("\n") ;

            if (nthreads != nthreads_max && 2 * nthreads > nthreads_max)
            {
                nthreads = nthreads_max ;
            }
            else
            {
                nthreads *= 2 ;
            }
        }
    }

    printf ("\n") ;
    LAGraph_FREE_ALL ;
    LAGraph_Finalize (NULL) ;
}
