//------------------------------------------------------------------------------
// msftest: test LAGraph_msf
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage: msftest can be used with both stdin or a file as its input.
// We assume by default that the matrix is symmetric. To override this,
// use the file-based input and pass 1 as the last argument.
//
// msftest < matrixmarketfile.mtx
// msftest matrixmarketfile.mtx
// msftest unsymmetric-matrixmarketfile.mtx 0
// msftest symmetric-matrixmarketfile.mtx 1

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING

#include <stdlib.h>  // for qsort
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&result);     \
    GrB_free (&A);          \
}

//****************************************************************************
GrB_Index *I, *J, *X;
int compare (const void * a, const void * b)
{
  return ( X[*(uint64_t*) a] - X[*(uint64_t*) b] );
}

//****************************************************************************
void check_solution (GrB_Matrix S, GrB_Matrix R)
{
    // check dimension
    GrB_Index svals, rvals, ns, nr;
    GrB_Matrix_nvals (&svals, S);
    GrB_Matrix_nvals (&rvals, R);
    GrB_Matrix_nrows (&ns, S);
    GrB_Matrix_nrows (&nr, R);
    if (ns != nr)
    {
        printf("wrong dimension\n");
        exit(0);
    }

    // check subset
    GrB_Index n = ns, cval, sum;
    GrB_Matrix U, C;
    GrB_Monoid Add;
    GrB_Monoid_new (&Add, GrB_PLUS_UINT64, (GrB_Index) 0);
    GrB_Matrix_new (&U, GrB_UINT64, n, n);
    GrB_Matrix_new (&C, GrB_UINT64, n, n);
    GrB_eWiseAdd (U, 0, 0, GrB_MIN_UINT64, S, R, 0);
    GrB_eWiseAdd (C, 0, 0, GxB_ISEQ_UINT64, U, S, 0);
    GrB_Matrix_nvals (&cval, C);
    GrB_reduce (&sum, 0, Add, C, 0);
    GrB_free (&U);
    GrB_free (&C);
    if (sum != cval || cval != svals)
    {
        printf("invalid set of edges\n");
        exit(0);
    }

    // check the spanning forest
    GrB_Index my_sol, answer = 0;
    GrB_reduce (&my_sol, 0, Add, R, 0);
    I = malloc (sizeof(GrB_Index) * svals);
    J = malloc (sizeof(GrB_Index) * svals);
    X = malloc (sizeof(GrB_Index) * svals);
    GrB_Matrix_extractTuples (I, J, X, &svals, S);
    GrB_Index *ind = malloc (sizeof(GrB_Index) * svals);
    GrB_Index *f = malloc (sizeof(GrB_Index) * n);
    for (GrB_Index i = 0; i < svals; i++)
        ind[i] = i;
    for (GrB_Index i = 0; i < n; i++)
        f[i] = i;
    qsort (ind, svals, sizeof(uint64_t), compare);
    for (GrB_Index i = 0; i < svals; i++)
    {
        GrB_Index x = I[ind[i]];
        GrB_Index y = J[ind[i]];
        bool comb = false;
        while (f[x] != f[y])
        {
            if (f[x] > f[y])
            {
                if (f[x] == x)
                {
                    comb = true;
                    f[x] = f[y];
                    break;
                }
                GrB_Index t = f[x];
                f[x] = f[y];
                x = t;
            }
            else
            {
                if (f[y] == y)
                {
                    comb = true;
                    f[y] = f[x];
                    break;
                }
                GrB_Index t = f[y];
                f[y] = f[x];
                y = t;
            }
        }
        if (comb)
            answer += X[ind[i]];
    }
    if (answer != my_sol)
    {
        printf("wrong answer!\n");
        printf("expected : %lu\n", answer);
        printf("actual   : %lu\n", my_sol);
        exit(0);
    }
    else
    {
        printf("correct (sum = %lu)\n", sum);
    }
    free(I); free(J); free(X); free(ind); free(f);
    GrB_free (&Add);
}

//****************************************************************************
int main (int argc, char **argv)
{
    GrB_Info info;
    GrB_Matrix A = NULL, S = NULL, result = NULL;
    GrB_Type A_type = NULL;
    GrB_init (GrB_NONBLOCKING);
    //LAGRAPH_OK (GxB_set (GxB_FORMAT, GxB_BY_ROW));

    FILE *f ;
    int symm = 0;
    if (argc == 1)
    {
        f = stdin ;
    }
    else
    {
        f = fopen (argv[1], "r");
        if (f == NULL)
        {
            printf ("unable to open file [%s]\n", argv[1]);
            return (GrB_INVALID_VALUE);
        }
        if (argc > 2)
            symm = atoi(argv[2]);
    }

    GrB_Index n;
    LAGRAPH_OK (LAGraph_MMRead (&A, &A_type, f, NULL));
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A));

    LAGRAPH_OK (GrB_Matrix_new (&S, GrB_UINT64, n, n));
    LAGRAPH_OK (GrB_eWiseAdd (S, 0, 0, GrB_MIN_UINT64, A, A, GrB_DESC_T1));

    #define NTRIALS 5
    int nthreads_max;
    int nthread_list [NTRIALS] = { 1, 4, 16, 20, 40 };

    LAGRAPH_OK (LAGraph_GetNumThreads(&nthreads_max, NULL));

    double tic[2], t;

    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        int nthreads = nthread_list [trial];
        if (nthreads > nthreads_max) break ;
        LAGraph_SetNumThreads (nthreads, NULL);
        printf("number of threads: %d\n", nthreads);

        LAGraph_Tic(tic, NULL);
        LAGRAPH_OK (LAGraph_msf (&result, S, true));
        LAGraph_Toc(&t, tic, NULL);
        check_solution (S, result);

        printf("Boruvka MSF: %f\n", t);
        printf("\n");
    }

    LAGraph_FREE_ALL ;
    LAGRAPH_OK (GrB_finalize ( ));
}
