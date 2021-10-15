//------------------------------------------------------------------------------
// LAGraph/Test/LCC/lcctest.c: test program for LAGraph_lcc
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M and Gabor Szarnyas, BME

// Usage: LCC_test can be used with both stdin or a file as its input.
// We assume by default that the matrix is symmetric. To override this,
// use the file-based input and pass 1 as the last argument.
//
// LCC_test < matrixmarketfile.mtx
// LCC_test matrixmarketfile.mtx
// LCC_test unsymmetric-matrixmarketfile.mtx 0
// LCC_test symmetric-matrixmarketfile.mtx 1

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL                            \
{                                                   \
    GrB_free (&C) ;                                 \
    GrB_free (&A) ;                                 \
    GrB_free (&M) ;                                 \
    GrB_free (&LCC) ;                               \
    GrB_free (&LCC1) ;                              \
    GrB_free (&LAGraph_ONE_FP64) ;                  \
}

//****************************************************************************

#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

void LAGraph_one_fp64
(
    double *z,
    const double *x       // ignored
)
{
    (*z) = 1 ;
}

//****************************************************************************
int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, C = NULL, M = NULL ;
    GrB_Vector LCC = NULL, LCC1 = NULL ;
    GrB_Type C_type = NULL, LCC_type = NULL;
    GrB_UnaryOp LAGraph_ONE_FP64 = NULL;

    LAGraph_Init (NULL) ;
    int nthreads_max;
    LAGraph_GetNumThreads (&nthreads_max, NULL) ;
    if (nthreads_max == 0) nthreads_max = 1 ;

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_Tic (tic, NULL) ;

    FILE *out = stdout ;

    FILE *f ;
    bool symmetric ;
    if (argc == 1)
    {
        f = stdin ;
        symmetric = false ;
    }
    else
    {
        f = fopen (argv[1], "r") ;
        if (f == NULL)
        {
            printf ("unable to open file [%s]\n", argv[1]) ;
            return (GrB_INVALID_VALUE) ;
        }
        if (argc == 2)
        {
            symmetric = false ;
        }
        else
        {
            symmetric = atoi(argv[2]) == 0;
        }
    }

    LAGRAPH_OK (LAGraph_MMRead (&C, &C_type, f, NULL)) ;
    GrB_Index n, ne ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, C)) ;
    LAGRAPH_OK (GrB_Matrix_nvals (&ne, C)) ;
    double t_read;
    LAGraph_Toc (&t_read, tic, NULL) ;
    fprintf (out, "\nread A time:     %14.6f sec\n", t_read) ;

    LAGraph_Tic (tic, NULL) ;

#if 1
#if LG_SUITESPARSE
    LAGraph_ONE_FP64   = GxB_ONE_FP64 ;
#else
    // create a new built-in unary operator using LAGraph_one_fp64
    LAGRAPH_OK (GrB_UnaryOp_new (&LAGraph_ONE_FP64,
                                 F_UNARY (LAGraph_one_fp64),
                                 GrB_FP64, GrB_FP64)) ;
#endif

    // A = spones (C), and typecast to FP64
    LAGRAPH_OK (GrB_Matrix_new (&A, GrB_FP64, n, n)) ;
    LAGRAPH_OK (GrB_apply (A, NULL, NULL, LAGraph_ONE_FP64, C, NULL)) ;
    GrB_free (&C) ;

    // M = diagonal mask matrix
    LAGRAPH_OK (GrB_Matrix_new (&M, GrB_BOOL, n, n)) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        // M(i,i) = true ;
        LAGRAPH_OK (GrB_Matrix_setElement (M, (bool) true, i, i)) ;
    }

    // remove self edges (via M)
    LAGRAPH_OK (GrB_assign (A, M, NULL, A, GrB_ALL, n, GrB_ALL, n,
                            GrB_DESC_RC));
    GrB_free (&M) ;

    LAGRAPH_OK (GrB_Matrix_nvals (&ne, A)) ;

    double t_process;
    LAGraph_Toc (&t_process, tic, NULL) ;
    fprintf (out, "process A time:  %14.6f sec\n", t_process) ;

#else
    A = C ;
    C = NULL ;
#endif

    LAGRAPH_OK (GrB_Matrix_nvals (&ne, A)) ;
#if LG_SUITESPARSE
    // GxB_fprint (A, GxB_SUMMARY, out) ;
#endif
    fprintf (out, "Matrix n: %0.16g, ne: %0.16g\n", (double) n, (double) ne) ;
    fflush (out) ;

    //--------------------------------------------------------------------------
    // compute LCC
    //--------------------------------------------------------------------------

    #define NTRIALS 5
    int nthread_list [NTRIALS] = { 1, 8, 16, 20, 40 } ;

    double t1 = -1 ;
    int nthreads_t1 = 0 ;
    // for (int nthreads = 1 ; nthreads <= nthreads_max ; )
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {

        int nthreads = nthread_list [trial] ;
        if (nthreads > nthreads_max) break ;
        LAGraph_SetNumThreads (nthreads, NULL) ;

        // ignore the sanitize time;  assume the user could have provided an
        // input graph that is already binary with no self-edges
        double timing [2] ;
        LAGRAPH_OK (LAGraph_lcc (&LCC, &LCC_type, A, symmetric, true, timing)) ;
        double t = timing [1] ;

        if (LCC1 == NULL)
        {
            LCC1 = LCC ;
            LCC = NULL ;
            t1 = t ;
            nthreads_t1 = nthreads ;
            // dump the result to lcc_results (required for comparing the
            // results with MATLAB)
#if 0
            FILE *results = fopen ("lcc_results", "w") ;
            for (GrB_Index i = 0 ; i < n ; i++)
            {
                double x = 0 ;
                LAGRAPH_OK (GrB_Vector_extractElement (&x, LCC1, i)) ;
                if (info == GrB_NO_VALUE) fprintf (results, " 0.\n") ;
                else fprintf (results, "%32.16g\n", x) ;
            }
            fclose (results) ;
#endif
        }
        else if (LCC1 != NULL)
        {
            bool ok ;
            LAGRAPH_OK (LAGraph_Vector_IsEqual_type (&ok, LCC, LCC1, LCC_type, NULL)) ;
            if (!ok) { fprintf (out, "error!\n") ; abort ( ) ; }
            GrB_free (&LCC) ; LCC = NULL;
        }

        fprintf (out, "nthreads: %3d sanitize %12.2f sec, LCC time: %10.2f"
            " sec, rate: %6.2f", nthreads, timing [0], t, 1e-6 * ne / t) ;
        if (nthreads != nthreads_t1 && t1 > 0)
        {
            fprintf (out, " speedup: %6.2f vs %d thread", t1 / t, nthreads_t1);
            if (nthreads_t1 != 1) fprintf (out, "s") ;
        }
        fprintf (out, "\n") ;
        fflush (out) ;
    }

    fprintf (out, "\n") ;
    LAGraph_FREE_ALL ;
    LAGraph_Finalize (NULL) ;
}
