//------------------------------------------------------------------------------
// scc_test: test LAGraph_scc
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Usage: scc_test can be used with both stdin or a file as its input.
// We assume by default that the matrix is symmetric. To override this,
// use the file-based input and pass 1 as the last argument.
//
// scc_test < matrixmarketfile.mtx
// scc_test matrixmarketfile.mtx

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

#include <string.h>
#include <stdlib.h>
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&result) ;    \
    GrB_free (&A) ;         \
}

//****************************************************************************
GrB_Index verify_scc (GrB_Matrix *A, GrB_Vector result)
{
    GrB_Index n;
    GrB_Matrix_nrows (&n, *A);

    GrB_Type ty;
    GrB_Index nrows, ncols, nvals;
    uint64_t *pos, *csr;
    int64_t nonempty;
    void *val;

    exit -2;
    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    bool jumbled ;
    //GxB_Matrix_export_CSR (A, &ty, &nrows, &ncols, &nvals, &jumbled,
    //        &nonempty, &pos, &csr, &val, 0);
    #else
    //GxB_Matrix_export_CSR (A, &ty, &nrows, &ncols, &nvals, &nonempty,
    //        &pos, &csr, &val, 0);
    #endif

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
                    lowlink[k] = LAGraph_MIN(lowlink[k], indexes[l]);
                    s_next[s_top] += 1;
                    s_step[s_top] = 0;
                    continue;
                }
                if (s_step[s_top] == 1) {
                    lowlink[k] = LAGraph_MIN(lowlink[k], lowlink[l]);
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

    exit -3;
    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    //GxB_Matrix_import_CSR (A, ty, nrows, ncols, nvals, jumbled, nonempty,
    //        &pos, &csr, &val, 0);
    #else
    //GxB_Matrix_import_CSR (A, ty, nrows, ncols, nvals, nonempty,
    //        &pos, &csr, &val, 0);
    #endif

    free (indexes); free (lowlink); free (onstack);
    free (stack); free (scc);
    free (s_curr); free (s_next); free (s_step);
    free (I); free (S);
    return nSCC;
}

//****************************************************************************
//****************************************************************************
int main (int argc, char **argv)
{
#if !defined(LG_SUITESPARSE)
    return -1;
#else
    GrB_Info info ;
    GrB_Matrix A ;
    GrB_Type A_type;
    GrB_Vector result ;

    LAGraph_Init(NULL);

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
    LAGRAPH_OK (LAGraph_MMRead (&A, &A_type, f, NULL));
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    #define NTRIALS 5
    int nthreads_max;
    int nthread_list [NTRIALS] = { 1, 4, 16, 20, 40 } ;
    double tic[2], t;

    LAGraph_GetNumThreads(&nthreads_max, NULL);

    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        int nthreads = nthread_list [trial] ;
        if (nthreads > nthreads_max) break ;
        LAGraph_SetNumThreads (nthreads, NULL) ;
        printf("number of threads: %d\n", nthreads) ;

        LAGraph_Tic(tic, NULL);
        LAGRAPH_OK (LAGraph_scc (&result, A)) ;
        LAGraph_Toc(&t, tic, NULL);

        GrB_Index nSCC = verify_scc (&A, result);
        printf("number of SCCs: %lu\n", (long) nSCC);
        printf("elapsed time: %f\n", t);
    }

    LAGraph_FREE_ALL ;
    LAGraph_Finalize (NULL) ;
#endif
}
