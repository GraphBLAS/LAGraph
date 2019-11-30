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

/**
 * Code is based on Borilvka's minimum spanning forest algorithm
 * Author: Yongzhe Zhang
 **/

#define LAGRAPH_FREE_ALL

#include "LAGraph.h"

static GrB_Info Matrix_removeElements_CSR (GrB_Matrix *A, GrB_Vector f)
{
    GrB_Info info;
    GrB_Index nrows, ncols, nvals;
    LAGRAPH_OK (GrB_Matrix_nrows (&nrows, *A));
    LAGRAPH_OK (GrB_Matrix_ncols (&ncols, *A));
    LAGRAPH_OK (GrB_Matrix_nvals (&nvals, *A));
    GrB_Index n = nrows;

    GrB_Index *csr, *pos;
    void *val;
    GrB_Type ty;
    int64_t nonempty = -1;
    // time consuming if A is not stored in CSR format
    GxB_Matrix_export_CSR (A, &ty, &nrows, &ncols, &nvals, &nonempty,
            &pos, &csr, &val, NULL);

    GrB_Index *I = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    GrB_Index *V = (GrB_Index*) malloc(sizeof(GrB_Index) * n);
    GrB_Vector_extractTuples(I, V, &n, f);
    int nthreads = LAGraph_get_nthreads();

    int64_t *range = (int64_t*) malloc (sizeof(int64_t) * (nthreads + 1));
    int64_t *count = (int64_t*) malloc (sizeof(int64_t) * (nthreads + 1));
    range[0] = 0;
    for (int i = 0; i < nthreads; i++)
        range[i + 1] = range[i] + (n + i) / nthreads;

    #pragma omp parallel for num_threads(nthreads)
    for (int id = 0; id < nthreads; id++) {
        int64_t ptr = pos[range[id]];
        for (int64_t v = range[id]; v < range[id + 1]; v++) {
            int64_t pv = V[v], start = pos[v];
            pos[v] = ptr;
            for (int64_t i = start; i < pos[v + 1]; i++) {
                int64_t u = csr[i];
                if (V[u] != V[v])
                    csr[ptr++] = u;
            }
        }
        count[id] = ptr - pos[range[id]];
    }

    int64_t offset = 0;
    for (int i = 0; i < nthreads; i++) {
        memcpy(csr + offset, csr + pos[range[i]], sizeof(int64_t) * count[i]);
        offset += count[i];
        count[i] = offset - count[i];
    }
    free(I); free(V);

    #pragma omp parallel for num_threads(nthreads)
    for (int id = 0; id < nthreads; id++) {
        int64_t ptr = pos[range[id]];
        for (int64_t v = range[id]; v < range[id + 1]; v++)
            pos[v] -= ptr - count[id];
    }
    pos[n] = offset;
    free(count); free(range);

    LAGRAPH_OK( GxB_Matrix_import_CSR (A, GrB_BOOL, n, n, offset, -1,
            &pos, &csr, &val, 0));
    return GrB_SUCCESS;
}

static GrB_Info Reduce_assign (GrB_Vector w, GrB_Vector mask, GrB_Vector src, GrB_Index *index, GrB_Index n)
{
    GrB_Index *mem = (GrB_Index*) malloc(sizeof(GrB_Index) * n * 4);
    GrB_Index *ind = mem, *sval = mem + n, *wval = sval + n, *mval = wval + n;
    GrB_Vector_extractTuples(ind, wval, &n, w);
    GrB_Vector_extractTuples(ind, sval, &n, src);
    GrB_Vector_extractTuples(ind, mval, &n, mask);
    for (GrB_Index i = 0; i < n; i++)
        if (mval[index[i]] && sval[i] < wval[index[i]])
            wval[index[i]] = sval[i];
    GrB_Vector_clear(w);
    GrB_Vector_build(w, ind, wval, n, GrB_PLUS_UINT64);
    free(mem);
    return GrB_SUCCESS;
}

GrB_Info LAGraph_cc_boruvka
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{
    GrB_Info info;
    GrB_Index n;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;

    GrB_Matrix S = NULL;
    if (sanitize)
    {
        GrB_Descriptor desc = NULL ;
        LAGRAPH_OK(GrB_Descriptor_new(&desc)) ;
        LAGRAPH_OK(GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN)) ;

        LAGRAPH_OK (GrB_Matrix_new (&S, GrB_BOOL, n, n)) ;
        LAGRAPH_OK (GrB_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, desc)) ;
        LAGRAPH_FREE(desc) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary and symmetric
        S = A ;
    }
    GrB_Vector f, p, m, i, e;
    LAGRAPH_OK (GrB_Vector_new (&f, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&p, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&m, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&i, GrB_UINT64, n));
    LAGRAPH_OK (GrB_Vector_new (&e, GrB_UINT64, n));
    // prepare
    GrB_Index *I = malloc(sizeof(GrB_Index) * n);
    GrB_Index *V = malloc(sizeof(GrB_Index) * n);
    for (GrB_Index i = 0; i < n; i++)
        I[i] = V[i] = i;
    LAGRAPH_OK (GrB_Vector_build (f, I, V, n, GrB_PLUS_UINT64));
    LAGRAPH_OK (GrB_assign (i, 0, 0, f, GrB_ALL, 0, 0));
    // semiring & monoid
    GrB_Monoid Min, Add;
    GrB_Semiring sel2ndMin;
    LAGRAPH_OK (GrB_Monoid_new (&Min, GrB_MIN_UINT64, n));
    LAGRAPH_OK (GrB_Semiring_new (&sel2ndMin, Min, GrB_SECOND_UINT64));
    LAGRAPH_OK (GrB_Monoid_new (&Add, GrB_PLUS_UINT64, (GrB_Index) 0));
    GrB_Index rem, diff;
    LAGRAPH_OK (GrB_Matrix_nvals (&rem, A));
    for (int iters = 1; rem > 0; iters++) {
        GrB_eWiseMult (e, 0, 0, GxB_ISEQ_UINT64, f, i, 0);
        GrB_assign (m, 0, 0, n, GrB_ALL, 0, 0); // m(:) = n
        GrB_mxv (m, 0, GrB_MIN_UINT64, sel2ndMin, S, f, 0);
        GrB_assign (p, 0, 0, n, GrB_ALL, 0, 0); // p(:) = n
        GrB_Vector_extractTuples (I, V, &n, f);
        Reduce_assign (p, e, m, V, n); // p[f[i]] = m[i] if e[f[i]]
        GrB_assign (m, 0, 0, n, GrB_ALL, 0, 0); // m(:) = n
        GrB_eWiseMult(e, 0, 0, GxB_ISNE_UINT64, m, p, 0); // e[i] = (p[i] != n)
        GrB_assign (f, e, 0, p, GrB_ALL, 0, 0); // f[i] = p[i] if (p[i] != n)
        GrB_Vector_extractTuples (I, V, &n, f);
        GrB_extract (p, 0, 0, f, V, n, 0); // p[i] = f[f[i]]
        GrB_eWiseMult (e, 0, 0, GxB_ISEQ_UINT64, p, i, 0);
        GrB_assign (p, 0, GrB_MIN_UINT64, f, GrB_ALL, 0, 0);
        GrB_assign (f, e, 0, p, GrB_ALL, 0, 0);
        do {
            GrB_Vector_extractTuples (I, V, &n, f);
            GrB_extract (p, 0, 0, f, V, n, 0);
            GrB_eWiseMult (e, 0, 0, GxB_ISNE_UINT64, f, p, 0);
            GrB_assign (f, 0, 0, p, GrB_ALL, 0, 0);
            GrB_reduce (&diff, 0, Add, e, 0);
        } while (diff != 0);
        Matrix_removeElements_CSR (&S, f);
        GrB_Matrix_nvals (&rem, S);
    }
    *result = f;

    free(I);
    free(V);
    GrB_free (&p);
    GrB_free (&i);
    GrB_free (&e);
    GrB_free (&m);
    GrB_free (&Add);
    GrB_free (&Min);
    GrB_free (&sel2ndMin);
    return GrB_SUCCESS;
}
