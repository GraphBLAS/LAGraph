//------------------------------------------------------------------------------
// scctest: test LAGraph_scc
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

// Usage: cctest can be used with both stdin or a file as its input.
// We assume by default that the matrix is symmetric. To override this,
// use the file-based input and pass 1 as the last argument.
//
// scctest < matrixmarketfile.mtx
// scctest matrixmarketfile.mtx

#include "LAGraph.h"
#include <string.h>
#include <stdlib.h>
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

#define MIN(a, b) ((a) < (b) ? (a) : (b))

GrB_Index verify_scc (GrB_Matrix *A, GrB_Vector result)
{
    GrB_Index n;
    GrB_Matrix_nrows (&n, *A);

    GrB_Type ty;
    GrB_Index nrows, ncols, nvals;
    uint64_t *pos, *csr;
    int64_t nonempty;
    void *val;
    GxB_Matrix_export_CSR (A, &ty, &nrows, &ncols, &nvals, &nonempty,
            &pos, &csr, &val, 0);

    int64_t *indexes = malloc (sizeof(int64_t) * n);
    int64_t *lowlink = malloc (sizeof(int64_t) * n);
    int64_t *onstack = malloc (sizeof(int64_t) * n);
    int64_t *stack = malloc (sizeof(int64_t) * n);
    int64_t *scc = malloc (sizeof(int64_t) * n);
    for (int64_t i = 0; i < n; i++)
        indexes[i] = -1;
    // stack simulation
    int64_t *s_curr = malloc (sizeof(int64_t) * n);
    int64_t *s_next = malloc (sizeof(int64_t) * n);
    int64_t *s_step = malloc (sizeof(int64_t) * n);
    int64_t s_top, index = 0, nSCC = 0, top = 0;

    for (int64_t i = 0; i < n; i++)
        if (indexes[i] == -1) {
            // push i
            s_top = 0;
            s_curr[0] = i;
            s_next[0] = pos[i];
            s_step[0] = 0;
            indexes[i] = index;
            lowlink[i] = index;
            index += 1;
            stack[top++] = i;
            onstack[i] = 1;
            // dfs
            while (s_top > -1) {
                int64_t k = s_curr[s_top];
                int64_t s = s_next[s_top];
                if (s == pos[k + 1]) {
                    // identify scc
                    if (indexes[k] == lowlink[k]) {
                        int64_t l;
                        do {
                            l = stack[--top];
                            scc[l] = k;
                            onstack[l] = 0;
                        } while (k != l);
                        nSCC++;
                    }
                    // pop
                    s_top--;
                    continue;
                }
                int64_t l = csr[s];
                if (s_step[s_top] == 0 && indexes[l] == -1) {
                    s_step[s_top] = 1;
                    // push l
                    ++s_top;
                    s_curr[s_top] = l;
                    s_next[s_top] = pos[l];
                    s_step[s_top] = 0;
                    indexes[l] = index;
                    lowlink[l] = index;
                    index += 1;
                    stack[top++] = l;
                    onstack[l] = true;
                    continue;
                }
                if (s_step[s_top] == 0 && onstack[l]) {
                    lowlink[k] = MIN(lowlink[k], indexes[l]);
                    s_next[s_top] += 1;
                    s_step[s_top] = 0;
                    continue;
                }
                if (s_step[s_top] == 1) {
                    lowlink[k] = MIN(lowlink[k], lowlink[l]);
                    s_next[s_top] += 1;
                    s_step[s_top] = 0;
                    continue;
                }
                // others 
                s_next[s_top] += 1;
                s_step[s_top] = 0;
                continue;
            }
        }

    uint64_t *I = malloc (sizeof(uint64_t) * n);
    uint64_t *S = malloc (sizeof(uint64_t) * n);
    GrB_Index len;
    GrB_Vector_nvals (&len, result);
    if (len != n) {
        printf("incorrect vector length ..\n");
        exit(0);
    }
    GrB_Vector_extractTuples (I, S, &len, result);
    int same = !memcmp(S, scc, sizeof(uint64_t) * n);
    if (!same) {
        printf("wrong answer!\n");
        exit(0);
    }

    GxB_Matrix_import_CSR (A, ty, nrows, ncols, nvals, nonempty,
            &pos, &csr, &val, 0);

    free (indexes); free (lowlink); free (onstack);
    free (stack); free (scc);
    free (s_curr); free (s_next); free (s_step);
    free (I); free (S);
    return nSCC;
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A ;
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
        printf("filename: %s\n", argv[1]);
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

    #define NTRIALS 5
    int nthreads_max;
    int nthread_list [NTRIALS] = { 1, 4, 16, 20, 40 } ;
    struct timeval t1, t2;

    LAGRAPH_OK (GxB_get (GxB_NTHREADS, &nthreads_max)) ;

    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        int nthreads = nthread_list [trial] ;
        if (nthreads > nthreads_max) break ;
        LAGraph_set_nthreads (nthreads) ;
        printf("number of threads: %d\n", nthreads) ;

        gettimeofday (&t1, 0) ;
        LAGRAPH_OK (LAGraph_scc (&result, A)) ;
        gettimeofday (&t2, 0) ;

        GrB_Index nSCC = verify_scc (&A, result);
        printf("number of SCCs: %lu\n", (long) nSCC);
        printf("elapsed time: %f\n", to_sec(t1, t2));
    }

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (GrB_finalize ( )) ;
}
