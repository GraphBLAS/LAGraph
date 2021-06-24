//------------------------------------------------------------------------------
// CC_test.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage: CC_test can be used with both stdin or a file as its input.
// We assume by default that the matrix is symmetric. To override this,
// use the file-based input and pass 1 as the last argument.
//
// CC_test < matrixmarketfile.mtx
// CC_test matrixmarketfile.mtx
// CC_test unsymmetric-matrixmarketfile.mtx 0
// CC_test symmetric-matrixmarketfile.mtx 1

//------------------------------------------------------------------------------

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&result) ;    \
    GrB_free (&A) ;         \
}

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

//****************************************************************************
GrB_Index countCC (GrB_Vector f, GrB_Index n)
{
    GrB_Index nCC = 0;
    GrB_Index *w_val = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    GrB_Vector_extractTuples(NULL, w_val, &n, f);
    for (GrB_Index i = 0; i < n; i++)
        if (w_val[i] == i)
            nCC += 1;
    free(w_val);
    return nCC;
}

//****************************************************************************
int main (int argc, char **argv)
{
    GrB_Info info ;

    GrB_Matrix A = NULL;
    GrB_Type   A_type = NULL;
    GrB_Vector result = NULL ;
    LAGRAPH_OK (LAGraph_Init (NULL)) ;
    //LAGRAPH_OK (GxB_set (GxB_BURBLE, false)) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max;
    LAGRAPH_OK (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
    if (Nthreads [1] == 0)
    {
        // create thread list automatically
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test: ") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    FILE *f ;

    GrB_Index n;

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;

    if (argc > 1)
    {
        // Usage:
        //      ./cc_test matrixfile.mtx
        //      ./cc_test matrixfile.grb

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }
        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            FILE *f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Binary file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_binread (&A, &A_type, f)) ;
            fclose (f) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            FILE *f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_MMRead(&A, &A_type, f, NULL));
            fclose (f) ;
        }

    }
    else
    {
        // Usage:  ./CC_test < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_MMRead(&A, &A_type, stdin, NULL));
    }

    ///////////
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    GrB_Index nvals, svals ;
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, A)) ;
    printf ("# of nodes: %lu  # of edges: %lu\n", n, nvals) ;

    double tic [2], t1, ttt, ttrial ;

    #define NTRIALS 16
    printf ("# of trials: %d\n\n", NTRIALS) ;

    bool sanitize = true ;

    GrB_Index nCC;
    for (int trial = 1 ; trial <= nt ; trial++)
    {
        int nthreads = Nthreads [trial] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_OK (LAGraph_SetNumThreads(nthreads, NULL)) ;

        t1 = 0 ;
        for (int k = 0 ; k < NTRIALS ; k++)
        {
            LAGraph_Tic (tic, NULL) ;
            LAGRAPH_OK (LAGraph_cc_boruvka (&result, A, sanitize)) ;
            LAGraph_Toc (&ttrial, tic, NULL) ;
            t1 += ttrial ;
            printf ("boruvka: trial: %2d time: %10.4f sec\n", k, ttrial) ;
            nCC = countCC (result, n) ;
            LAGRAPH_OK (GrB_free (&result)) ;
        }
        ttt = t1 / NTRIALS ;
        printf("Boruvka: threads: %2d time: %10.4f  # of CC: %lu\n\n",
               nthreads, ttt, nCC) ;
        fprintf (stderr, "Avg: CC (boruvka) %3d: %10.3f sec: %s\n",
                 nthreads, ttt, matrix_name) ;

        t1 = 0 ;
        for (int k = 0 ; k < NTRIALS ; k++)
        {
            LAGraph_Tic (tic, NULL) ;
            LAGRAPH_OK (LAGraph_cc_lacc (&result, A, sanitize)) ;
            LAGraph_Toc (&ttrial, tic, NULL) ;
            t1 += ttrial ;
            printf ("lacc: trial: %2d time: %10.4f sec\n", k, ttrial) ;
            nCC = countCC (result, n) ;
            LAGRAPH_OK (GrB_free (&result)) ;
        }
        ttt = t1 / NTRIALS ;
        printf("LACC: threads: %2d time: %10.4f  # of CC: %lu\n\n",
               nthreads, ttt, nCC) ;
        fprintf (stderr,
                 "Avg: CC (lacc) %3d: %10.3f sec: %s\n",
                 nthreads, ttt, matrix_name) ;
    }

    printf ("\n") ;
    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (LAGraph_Finalize (NULL)) ;
}
