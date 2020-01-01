//------------------------------------------------------------------------------
// cctest: test LAGraph_cc_*.c
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

// Contributed by Tim Davis, Texas A&M

// Usage: cctest can be used with both stdin or a file as its input.
// We assume by default that the matrix is symmetric. To override this,
// use the file-based input and pass 1 as the last argument.
//
// cctest < matrixmarketfile.mtx
// cctest matrixmarketfile.mtx
// cctest unsymmetric-matrixmarketfile.mtx 0
// cctest symmetric-matrixmarketfile.mtx 1

#include "LAGraph.h"
#include <sys/time.h>

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&result) ;    \
    GrB_free (&A) ;         \
}

double to_sec(struct timeval t1, struct timeval t2)
{
    return
        (t2.tv_sec  - t1.tv_sec ) + 
        (t2.tv_usec - t1.tv_usec) * 1e-6;
}

GrB_Index countCC (GrB_Vector f, GrB_Index n)
{
    GrB_Index nCC = 0;
    GrB_Index *w_ind = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    GrB_Index *w_val = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    GrB_Vector_extractTuples(w_ind, w_val, &n, f);
    for (GrB_Index i = 0; i < n; i++)
        if (w_val[i] == i)
            nCC += 1;
    free(w_ind);
    free(w_val);
    return nCC;
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A, S ;
    GrB_Vector result ;
    GrB_init (GrB_NONBLOCKING) ;
    LAGRAPH_OK (GxB_set (GxB_FORMAT, GxB_BY_ROW)) ;

    FILE *f ;
    int symm = 0;
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
        if (argc > 2)
            symm = atoi(argv[2]);
    }

    GrB_Index n;
    LAGRAPH_OK (LAGraph_mmread (&A, f)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    GrB_Descriptor desc = 0 ;
    LAGRAPH_OK (GrB_Descriptor_new(&desc)) ;
    LAGRAPH_OK (GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN)) ;

    LAGRAPH_OK (GrB_Matrix_new (&S, GrB_BOOL, n, n)) ;
    LAGRAPH_OK (GrB_eWiseAdd (S, 0, 0, GrB_LOR, A, A, desc)) ;
    LAGRAPH_FREE (desc) ;

    #define NTRIALS 5
    int nthreads_max;
    int nthread_list [NTRIALS] = { 1, 4, 16, 20, 40 } ;
    struct timeval t1, t2;

    LAGRAPH_OK (GxB_get (GxB_NTHREADS, &nthreads_max)) ;

    GrB_Index nCC;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        int nthreads = nthread_list [trial] ;
        if (nthreads > nthreads_max) break ;
        LAGraph_set_nthreads (nthreads) ;
        printf("number of threads: %d\n", nthreads) ;

        gettimeofday (&t1, 0) ;
        LAGRAPH_OK (LAGraph_cc_fastsv (&result, A, true)) ;
        gettimeofday (&t2, 0) ;

        nCC = countCC (result, n) ;
        printf("number of CCs: %lu\n", nCC) ;
        printf("FastSV: %f\n", to_sec (t1, t2)) ;

        gettimeofday (&t1, 0) ;
        LAGRAPH_OK (LAGraph_cc_boruvka (&result, A, true)) ;
        gettimeofday (&t2, 0) ;

        nCC = countCC (result, n) ;
        printf("number of CCs: %lu\n", nCC) ;
        printf("Boruvka: %f\n", to_sec (t1, t2)) ;
        printf("\n");
    }

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (GrB_finalize ( )) ;
}

