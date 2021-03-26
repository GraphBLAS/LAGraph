//------------------------------------------------------------------------------
// LAGraph/Test/CDLP/cdlptest.c: test program for LAGraph_cdlp
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

// Contributed by Tim Davis, Texas A&M

// Usage:  cdlptest < matrixmarketfile.mtx
//         cdlptest matrixmarketfile.mtx

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#define LAGRAPH_FREE_ALL                            \
{                                                   \
    GrB_free (&C) ;                                 \
    GrB_free (&A) ;                                 \
    GrB_free (&M) ;                                 \
}

void Print_Label_Matrix(GrB_Matrix m)
{
    GrB_Index row_vals;
    GrB_Matrix_nrows(&row_vals, m);

    printf("Label vec:\n");

    uint64_t value;

    printf(" ");
    for (GrB_Index i = 0; i < row_vals; i++)
    {
        printf(" %ld", i + 1);
    }
    printf("\n");

    printf("[");
    for (GrB_Index i = 0; i < row_vals; i++)
    {
        GrB_Matrix_extractElement(&value, m, i, i);
        printf(" %ld", value);
    }
    printf(" ]\n");
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, C = NULL, M = NULL ;
    GrB_Vector CDLP = NULL, CDLP1 = NULL ;

    LAGraph_init ( ) ;
    int nthreads_max = LAGraph_get_nthreads ( ) ;
    if (nthreads_max == 0) nthreads_max = 1 ;

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;

    FILE *out = stdout ;

    FILE *f ;
    bool symmetric ;
    int itermax ;
    if (argc < 4)
    {
        printf ("Usage: matrixmarketfile.mtx is_symmetric max_number_of_iterations" ) ;
        return (GrB_INVALID_VALUE) ;
    }
    else
    {
        f = fopen (argv[1], "r") ;
        if (f == NULL)
        {
            printf ("unable to open file [%s]\n", argv[1]) ;
            return (GrB_INVALID_VALUE) ;
        }
        symmetric = atoi(argv[2]) ;
        itermax = atoi(argv[3]) ;
    }

    LAGRAPH_OK (LAGraph_mmread (&C, f)) ;
    GrB_Index n, ne ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, C)) ;
    LAGRAPH_OK (GrB_Matrix_nvals (&ne, C)) ;
    double t_read = LAGraph_toc (tic) ;
    fprintf (out, "\nread A time:     %14.6f sec\n", t_read) ;

    LAGraph_tic (tic) ;

#if 0
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
        LAGraph_desc_oocr)) ;
    GrB_free (&M) ;

    LAGRAPH_OK (GrB_Matrix_nvals (&ne, A)) ;

    // double t_process = LAGraph_toc (tic) ;
    fprintf (out, "process A time:  %14.6f sec\n", t_process) ;

#else
    A = C ;
    C = NULL ;
#endif

    LAGRAPH_OK (GrB_Matrix_nvals (&ne, A)) ;
    #ifdef GxB_SUITESPARSE_GRAPHBLAS
    // GxB_fprint (A, GxB_SUMMARY, out) ;
    #endif
    fprintf (out, "Matrix n: %0.16g, ne: %0.16g\n", (double) n, (double) ne) ;
    fflush (out) ;

    //--------------------------------------------------------------------------
    // compute CDLP
    //--------------------------------------------------------------------------

    #define NTRIALS 4
    int nthread_list [NTRIALS] = { 1, 2, 4, 8 };

    double t1 = -1 ;
    int nthreads_t1 = 0 ;
    // for (int nthreads = 1 ; nthreads <= nthreads_max ; )
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {

        int nthreads = nthread_list [trial] ;
        if (nthreads > nthreads_max) break ;
        LAGraph_set_nthreads (nthreads) ;

        // ignore the sanitize time;  assume the user could have provided an
        // input graph that is already binary with no self-edges
        double timing [2] ;
        LAGRAPH_OK (LAGraph_cdlp(&CDLP, A, symmetric, true, itermax, timing)) ;
        double t = timing [1] ;
//        GxB_print(CDLP, GxB_COMPLETE);

        if (CDLP1 == NULL)
        {
            CDLP1 = CDLP ;
            CDLP = NULL ;
            t1 = t ;
            nthreads_t1 = nthreads ;
        }
        else if (CDLP1 != NULL)
        {
            bool ok ;
            LAGRAPH_OK (LAGraph_Vector_isequal (&ok, CDLP, CDLP1, GrB_EQ_FP64)) ;
            if (!ok) { fprintf (out, "error!\n") ; abort ( ) ; }
            GrB_free (&CDLP) ;
        }

        fprintf (out, "nthreads: %3d sanitize %12.2f sec, CDLP time: %10.2f"
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
    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}

