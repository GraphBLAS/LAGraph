//------------------------------------------------------------------------------
// LAGraph_cdlp: community detection using label propagation
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

// LAGraph_cdlp: Contributed by Gabor Szarnyas and Balint Hegyi,
// Budapest University of Technology and Economics
// (with accented characters: G\'{a}bor Sz\'{a}rnyas and B\'{a}lint Hegyi,
// using LaTeX syntax). https://inf.mit.bme.hu/en/members/szarnyasg .

// ## Background
//
// This function was originally written for the LDBC Graphalytics benchmark.
//
// The community detection using label propagation (CDLP) algorithm is
// defined both for directed and undirected graphs.
//
// The definition implemented here is described in the following document:
// https://ldbc.github.io/ldbc_graphalytics_docs/graphalytics_spec.pdf
//
// The algorithm is based on the one given in the following paper:
//
// Usha Raghavan, Réka Albert, and Soundar Kumara. "Near linear time algorithm
// to detect community structures in large-scale networks". In: Physical
// Review E 76.3 (2007), p. 036106, https://arxiv.org/abs/0709.2938
//
// The key idea of the algorithm is that each vertex is assigned the label
// that is most frequent among its neighbors. To allow reproducible
// experiments, the algorithm is modified to guarantee deterministic behavior:
// it always picks the smallest label in case of a tie:
//
// min ( argmax_{l} (#neighbors with label l) )
//
// In other words, we need to compute the *minimum mode value* (minmode) for
// the labels among the neighbors.
//
// For directed graphs, a label on a neighbor that is connected through both
// an outgoing and on an incoming edge counts twice:
//
// min ( argmax_{l} (#incoming neighbors with l + #outgoing neighbors with l) )
//
// ## Example (undirected)
//
// For an example, let's assume an undirected graph where vertex 1 has four
// neighbors {2, 3, 4, 5}, and the current labels in the graph are
// L = [3, 5, 4, 5, 4].
//
// In this example, the distribution of labels among the neighbors of vertex 1
// is {4 => 2, 5 => 2}, therefore, the minimum mode value is 4.
//
// Next, we capture this operation using GraphBLAS operations and
// data structures. Notice that the neighbors of vertex 1 are encoded
// as a sparse vector in the adjacency matrix:
//
// A = | 0 1 1 1 1 |
//     | 1 . . .   |
//     | 1 .       |
//     | 1 .       |
//     | 1         |
//
// To allow propagating the labels along edges, we use a diagonal matrix
// with the elements of the diagonal set to the values of L:
//
// diag(L) = | 3 0 0 0 0 |
//           | 0 5 0 0 0 |
//           | 0 0 4 0 0 |
//           | 0 0 0 5 0 |
//           | 0 0 0 0 4 |
//
// If we multiply adjacency matrix with diag(L), we get a matrix
// containing the labels of the neighbor nodes.
//
// In the example, this gives:
//
// AL = A*diag(L) = | 0 5 4 5 4 |
//                  | . . .     |
//
// ## Selecting the minimum mode value
//
// Next, we need to compute the minimum mode value for each row. As it is
// difficult to capture this operation as a monoid, we use a sort operation
// on each row. In the undirected case, we extract tuples <I, _, X> from the
// matrix, then use <I, X> for sorting. In the directed case, we extract
// tuples <I1, _, X1> and <I2, _, X2>, then use <I1+I2, X1+X2>,
// where '+' denotes concatenation. Column indices (J) are not used.
//
// The resulting two-tuples are sorted using a parallel merge sort.
// Finally, we use the sorted arrays compute the minimum mode value for each
// row.
//
// ## Fixed point
//
// At the end of each iteration, we check whether L[i-1] == L[i] and
// terminate if we reached a fixed point.
//
// ## Further optimizations
//
// A possible optimization is that the first iteration is rather trivial:
//
// * in the undirected case, each vertex gets the minimal initial label (=id)
//   of its neighbors.
// * in the directed case, each vertex gets the minimal initial label (=id)
//   of its neighbors which are doubly-linked (on an incoming and on an
//   outgoing edge). In the absence of such a neighbor, it picks the minimal
//   label of its neighbors (connected through either an incoming or through
//   an outgoing edge).

#include "LAGraph_internal.h"
#include "GB_msort_2.h"

#define LAGRAPH_FREE_ALL                                                       \
{                                                                              \
    GrB_free (&L) ;                                                            \
    if (sanitize) GrB_free (&S) ;                                              \
    GrB_free (&S) ;                                                            \
    GrB_free (&AL_in) ;                                                        \
    GrB_free (&AL_out) ;                                                       \
    GrB_free (&desc_in) ;                                                      \
    GrB_free (&desc_out) ;                                                     \
}

GrB_Info LAGraph_cdlp
(
    GrB_Vector *CDLP_handle, // output vector
    const GrB_Matrix A,      // input matrix
    bool symmetric,          // denote whether the matrix is symmetric
    bool sanitize,           // if true, ensure A is binary
    int itermax,             // max number of iterations,
    double *t                // t [0] = sanitize time, t [1] = cdlp time,
                             // in seconds
)
{
    GrB_Info info;

    // Diagonal label matrix
    GrB_Matrix L = NULL;
    GrB_Matrix L_prev = NULL;
    // Source adjacency matrix
    GrB_Matrix S = NULL;
    // S*L matrix
    GrB_Matrix AL_in = NULL;
    GrB_Matrix AL_out = NULL;
    // Result CDLP vector
    GrB_Vector CDLP = NULL;

    GrB_Descriptor desc_in = NULL;
    GrB_Descriptor desc_out = NULL;

    // Arrays holding extracted tuples
    GrB_Index *I = NULL;
    GrB_Index *X = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (CDLP_handle == NULL)
    {
        return GrB_NULL_POINTER;
    }

    //--------------------------------------------------------------------------
    // ensure input is binary and has no self-edges
    //--------------------------------------------------------------------------

    double tic [2];
    t [0] = 0;         // sanitize time
    t [1] = 0;         // CDLP time

    if (sanitize)
    {
        LAGraph_tic (tic) ;

        // S = binary pattern of A
        LAGRAPH_OK (LAGraph_pattern(&S, A, GrB_UINT64))
        // Remove all self edges
        LAGRAPH_OK (LAGraph_prune_diag(S))

        t [0] = LAGraph_toc (tic) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary with no self edges.
        // Results are undefined if this condition does not hold.
        S = A;
    }

    LAGraph_tic (tic) ;

    GxB_Format_Value A_format = -1;
    LAGRAPH_OK (GxB_get(A, GxB_FORMAT, &A_format))
    if (A_format != GxB_BY_ROW)
    {
        LAGRAPH_ERROR(
            "CDLP algorithm only works on matrices stored by row (CSR)",
            GrB_INVALID_OBJECT
        )
    }

    LAGRAPH_OK(GrB_Descriptor_new(&desc_in))
    LAGRAPH_OK(GrB_Descriptor_set(desc_in, GrB_OUTP, GrB_REPLACE))

    LAGRAPH_OK(GrB_Descriptor_new(&desc_out))
    LAGRAPH_OK(GrB_Descriptor_set(desc_out, GrB_INP0, GrB_TRAN))
    LAGRAPH_OK(GrB_Descriptor_set(desc_out, GrB_OUTP, GrB_REPLACE))


    // n = size of A (# of nodes in the graph)
    GrB_Index n;
    LAGRAPH_OK (GrB_Matrix_nrows(&n, A))

    // nz = # of non-zero elements in the matrix
    // nnz = # of non-zero elements used in the computations
    //   (twice as many for directed graphs)
    GrB_Index nz, nnz;
    LAGRAPH_OK (GrB_Matrix_nvals(&nz, A))
    if (!symmetric)
    {
        nnz = 2 * nz;
    }
    else
    {
        nnz = nz;
    }

    // Initialize L with diagonal elements 1..n
    LAGRAPH_OK(GrB_Matrix_new(&L, GrB_UINT64, n, n))
    for (GrB_Index i = 0; i < n; i++)
    {
        LAGRAPH_OK(GrB_Matrix_setElement(L, i + 1, i, i))
    }
    // Initialize matrix for storing previous labels
    LAGRAPH_OK(GrB_Matrix_new(&L_prev, GrB_UINT64, n, n))

    LAGRAPH_OK(GrB_Matrix_new(&AL_in, GrB_UINT64, n, n))
    if (!symmetric)
    {
        LAGRAPH_OK(GrB_Matrix_new(&AL_out, GrB_UINT64, n, n))
    }

    // Initialize data structures for extraction from 'AL_in' and (for directed graphs) 'AL_out'
    I = LAGraph_malloc(nnz, sizeof(GrB_Index));
    X = LAGraph_malloc(nnz, sizeof(GrB_Index));

    uint64_t* workspace1 = LAGraph_malloc(nnz, sizeof(GrB_Index));
    uint64_t* workspace2 = LAGraph_malloc(nnz, sizeof(GrB_Index));

    const int nthreads = LAGraph_get_nthreads();
    for (int iteration = 0; iteration < itermax; iteration++)
    {
        // AL_in = A * L
        LAGRAPH_OK(GrB_mxm(AL_in, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, S, L, desc_in))
        LAGRAPH_OK(GrB_Matrix_extractTuples_UINT64(I, GrB_NULL, X, &nz, AL_in))

        if (!symmetric)
        {
            // AL_out = A' * L
            LAGRAPH_OK(GrB_mxm(AL_out, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, S, L, desc_out))
            LAGRAPH_OK(GrB_Matrix_extractTuples_UINT64(&I[nz], GrB_NULL, &X[nz], &nz, AL_out))
        }

        GB_msort_2(I, X, workspace1, workspace2, nnz, nthreads);

        // save current labels for comparison by swapping L and L_prev
        GrB_Matrix L_swap = L;
        L = L_prev;
        L_prev = L_swap;

        GrB_Index mode_value = -1;
        GrB_Index mode_length = 0;
        GrB_Index run_length = 1;

        // I[k] is the current row index
        // X[k] is the current value
        // we iterate in range 1..nnz and use the last index (nnz) to process the last row of the matrix
        for (GrB_Index k = 1; k <= nnz; k++)
        {
            // check if we have a reason to recompute the mode value
            if (k == nnz        // we surpassed the last element
             || I[k-1] != I[k]  // the run value has changed
             || X[k-1] != X[k]) // the row index has changed
            {
                if (run_length > mode_length)
                {
                    mode_value = X[k-1];
                    mode_length = run_length;
                }
                run_length = 0;
            }
            run_length++;

            // check if we passed a row
            if (k == nnz        // we reached the last element
             || I[k-1] != I[k]) // the row index has changed
            {
                GrB_Matrix_setElement(L, mode_value, I[k-1], I[k-1]);
                mode_length = 0;
            }
        }

        if (L_prev == L)
        {
            bool same = true;
            for (GrB_Index i = 0; i < n; i++)
            {
                uint64_t x_previous;
                GrB_Matrix_extractElement(&x_previous, L_prev, i, i);

                uint64_t x;
                GrB_Matrix_extractElement(&x, L, i, i);

                if (x != x_previous) same = false;
            }
            if (same)
            {
                break;
            }
        }
    }

    //--------------------------------------------------------------------------
    // extract final labels to the result vector
    //--------------------------------------------------------------------------

    LAGRAPH_OK (GrB_Vector_new(&CDLP, GrB_UINT64, n))
    for (GrB_Index i = 0; i < n; i++)
    {
        uint64_t x;
        LAGRAPH_OK(GrB_Matrix_extractElement(&x, L, i, i))
        LAGRAPH_OK(GrB_Vector_setElement(CDLP, x, i))
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*CDLP_handle) = CDLP;
    CDLP = NULL;            // set to NULL so LAGRAPH_FREE_ALL doesn't free it
    LAGRAPH_FREE_ALL

    t [1] = LAGraph_toc (tic) ;

    return (GrB_SUCCESS);
}
