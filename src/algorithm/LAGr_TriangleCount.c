//------------------------------------------------------------------------------
// LAGr_TriangleCount: Triangle counting using various methods
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Count the number of triangles in a graph,

// This is an Advanced algorithm (G->ndiag, G->rowdegree,
// G->structure_is_symmetric are required).

// Given a symmetric graph A with no-self edges, LAGr_TriangleCount counts the
// number of triangles in the graph.  A triangle is a clique of size three,
// that is, 3 nodes that are all pairwise connected.

// One of 6 methods are used, defined below where L and U are the strictly
// lower and strictly upper triangular parts of the symmetrix matrix A,
// respectively.  Each method computes the same result, ntri:

//  0:  default:    use the default method (currently methood 5)
//  1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia:     ntri = sum (sum ((L * L) .* L))
//  4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
//  5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  SandiaDot2: ntri = sum (sum ((U * L') .* U)).  Note that U=L'.

// A is a square symmetric matrix, of any type.  Its values are ignored.
// Results are undefined for methods 1 and 2 if self-edges exist in A.  Results
// are undefined for all methods if A is unsymmetric.

// The Sandia* methods all tend to be faster than the Burkhardt or Cohen
// methods.  For the largest graphs, SandiaDot tends to be fastest, except for
// the GAP-urand matrix, where the saxpy-based Sandia method (L*L.*L) is
// fastest.  For many small graphs, the saxpy-based Sandia and Sandia2 methods
// are often faster that the dot-product-based methods.

// Reference (for the "Sandia*" methods): Wolf, Deveci, Berry, Hammond,
// Rajamanickam, 'Fast linear algebra- based triangle counting with
// KokkosKernels', IEEE HPEC'17, https://dx.doi.org/10.1109/HPEC.2017.8091043,

#define LG_FREE_ALL             \
{                               \
    GrB_free (L) ;              \
    GrB_free (U) ;              \
}

#include "LG_internal.h"

//------------------------------------------------------------------------------
// tricount_prep: construct L and U for LAGr_TriangleCount
//------------------------------------------------------------------------------

static int tricount_prep
(
    GrB_Matrix *L,      // if present, compute L = tril (A,-1)
    GrB_Matrix *U,      // if present, compute U = triu (A, 1)
    GrB_Matrix A,       // input matrix
    char *msg
)
{
    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    if (L != NULL)
    {
        // L = tril (A,-1)
        GRB_TRY (GrB_Matrix_new (L, GrB_BOOL, n, n)) ;
        GRB_TRY (GrB_select (*L, NULL, NULL, GrB_TRIL, A, (int64_t) (-1),
            NULL)) ;
    }

    if (U != NULL)
    {
        // U = triu (A,1)
        GRB_TRY (GrB_Matrix_new (U, GrB_BOOL, n, n)) ;
        GRB_TRY (GrB_select (*U, NULL, NULL, GrB_TRIU, A, (int64_t) 1, NULL)) ;
    }
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
{                                           \
    GrB_free (&C) ;                         \
    GrB_free (&L) ;                         \
    GrB_free (&T) ;                         \
    GrB_free (&U) ;                         \
    LAGraph_Free ((void **) &P, NULL) ;     \
}

int LAGr_TriangleCount
(
    // output:
    uint64_t       *ntriangles,
    // input:
    const LAGraph_Graph G,
    LAGraph_TriangleCount_Method    method,
    LAGraph_TriangleCount_Presort *presort,
    char           *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL ;
    int64_t *P = NULL ;
    LG_ASSERT_MSG (
    method == LAGraph_TriangleCount_Default ||   // 0: use default method
    method == LAGraph_TriangleCount_Burkhardt || // 1: sum (sum ((A^2) .* A))/6
    method == LAGraph_TriangleCount_Cohen ||     // 2: sum (sum ((L * U) .*A))/2
    method == LAGraph_TriangleCount_Sandia ||    // 3: sum (sum ((L * L) .* L))
    method == LAGraph_TriangleCount_Sandia2 ||   // 4: sum (sum ((U * U) .* U))
    method == LAGraph_TriangleCount_SandiaDot || // 5: sum (sum ((L * U') .* L))
    method == LAGraph_TriangleCount_SandiaDot2,  // 6: sum (sum ((U * L') .* U))
    GrB_INVALID_VALUE, "method is invalid") ;
    if (presort != NULL)
    {
        LG_ASSERT_MSG (
        (*presort) == LAGraph_TriangleCount_NoSort ||
        (*presort) == LAGraph_TriangleCount_Ascending ||
        (*presort) == LAGraph_TriangleCount_Descending ||
        (*presort) == LAGraph_TriangleCount_AutoSort,
        GrB_INVALID_VALUE, "presort is invalid") ;
    }
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (ntriangles != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (G->ndiag == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED) ;

    if (method == LAGraph_TriangleCount_Default)
    {
        // 0: use default method (5): SandiaDot: sum (sum ((L * U') .* L))
        method = LAGraph_TriangleCount_SandiaDot ;
    }

    LG_ASSERT_MSG ((G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->structure_is_symmetric == LAGraph_TRUE)),
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be known to be symmetric") ;

    // the Sandia* methods can benefit from the presort
    bool method_can_use_presort =
    method == LAGraph_TriangleCount_Sandia ||    // 3: sum (sum ((L * L) .* L))
    method == LAGraph_TriangleCount_Sandia2 ||   // 4: sum (sum ((U * U) .* U))
    method == LAGraph_TriangleCount_SandiaDot || // 5: sum (sum ((L * U') .* L))
    method == LAGraph_TriangleCount_SandiaDot2 ; // 6: sum (sum ((U * L') .* U))

    GrB_Matrix A = G->A ;
    GrB_Vector Degree = G->rowdegree ;
    bool auto_sort = (presort != NULL)
        && ((*presort) == LAGraph_TriangleCount_AutoSort) ;
    if (auto_sort && method_can_use_presort)
    {
        LG_ASSERT_MSG (Degree != NULL,
            LAGRAPH_PROPERTY_MISSING, "G->rowdegree is required") ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;
    GrB_Semiring semiring = LAGraph_plus_one_int64 ;
    GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;

    //--------------------------------------------------------------------------
    // heuristic sort rule
    //--------------------------------------------------------------------------

    if (auto_sort)
    {
        // auto selection of sorting method
        (*presort) = LAGraph_TriangleCount_NoSort ; // default is not to sort

        if (method_can_use_presort)
        {
            // This rule is very similar to Scott Beamer's rule in the GAP TC
            // benchmark, except that it is extended to handle the ascending
            // sort needed by methods 3 and 5.  It also uses a stricter rule,
            // since the performance of triangle counting in SuiteSparse:
            // GraphBLAS is less sensitive to the sorting as compared to the
            // GAP algorithm.  This is because the dot products in SuiteSparse:
            // GraphBLAS use binary search if one vector is very sparse
            // compared to the other.  As a result, SuiteSparse:GraphBLAS needs
            // the sort for fewer matrices, as compared to the GAP algorithm.

            // With this rule, the GAP-kron and GAP-twitter matrices are
            // sorted, and the others remain unsorted.  With the rule in the
            // GAP tc.cc benchmark, GAP-kron and GAP-twitter are sorted, and so
            // is GAP-web, but GAP-web is not sorted here.

            #define NSAMPLES 1000
            GrB_Index nvals ;
            GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
            if (n > NSAMPLES && ((double) nvals / ((double) n)) >= 10)
            {
                // estimate the mean and median degrees
                double mean, median ;
                LG_TRY (LAGraph_SampleDegree (&mean, &median,
                    G, true, NSAMPLES, n, msg)) ;
                // sort if the average degree is very high vs the median
                if (mean > 4 * median)
                {
                    switch (method)
                    {
                        case LAGraph_TriangleCount_Sandia:
                            // 3:sum (sum ((L * L) .* L))
                            (*presort) = LAGraph_TriangleCount_Ascending  ;
                            break ;
                        case LAGraph_TriangleCount_Sandia2:
                            // 4: sum (sum ((U * U) .* U))
                            (*presort) = LAGraph_TriangleCount_Descending ;
                            break ;
                        default:
                        case LAGraph_TriangleCount_SandiaDot:
                            // 5: sum (sum ((L * U') .* L))
                            (*presort) = LAGraph_TriangleCount_Ascending  ;
                            break ;
                        case LAGraph_TriangleCount_SandiaDot2:
                            // 6: sum (sum ((U * L') .* U))
                            (*presort) = LAGraph_TriangleCount_Descending ;
                            break ;
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // sort the input matrix, if requested
    //--------------------------------------------------------------------------

    if (presort != NULL && (*presort) != LAGraph_TriangleCount_NoSort)
    {
        // P = permutation that sorts the rows by their degree
        LG_TRY (LAGraph_SortByDegree (&P, G, true,
            (*presort) == LAGraph_TriangleCount_Ascending, msg)) ;

        // T = A (P,P) and typecast to boolean
        GRB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;
        GRB_TRY (GrB_extract (T, NULL, NULL, A, (GrB_Index *) P, n,
            (GrB_Index *) P, n, NULL)) ;
        A = T ;

        // free workspace
        LG_TRY (LAGraph_Free ((void **) &P, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // count triangles
    //--------------------------------------------------------------------------

    int64_t ntri ;

    switch (method)
    {

        case LAGraph_TriangleCount_Burkhardt:  // 1: sum (sum ((A^2) .* A)) / 6

            GRB_TRY (GrB_mxm (C, A, NULL, semiring, A, A, GrB_DESC_S)) ;
            GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            ntri /= 6 ;
            break ;

        case LAGraph_TriangleCount_Cohen: // 2: sum (sum ((L * U) .* A)) / 2

            LG_TRY (tricount_prep (&L, &U, A, msg)) ;
            GRB_TRY (GrB_mxm (C, A, NULL, semiring, L, U, GrB_DESC_S)) ;
            GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            ntri /= 2 ;
            break ;

        case LAGraph_TriangleCount_Sandia: // 3: sum (sum ((L * L) .* L))

            // using the masked saxpy3 method
            LG_TRY (tricount_prep (&L, NULL, A, msg)) ;
            GRB_TRY (GrB_mxm (C, L, NULL, semiring, L, L, GrB_DESC_S)) ;
            GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case LAGraph_TriangleCount_Sandia2: // 4: sum (sum ((U * U) .* U))

            // using the masked saxpy3 method
            LG_TRY (tricount_prep (NULL, &U, A, msg)) ;
            GRB_TRY (GrB_mxm (C, U, NULL, semiring, U, U, GrB_DESC_S)) ;
            GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        default:
        case LAGraph_TriangleCount_SandiaDot: // 5: sum (sum ((L * U') .* L))

            // This tends to be the fastest method for most large matrices, but
            // the SandiaDot2 method is also very fast.

            // using the masked dot product
            LG_TRY (tricount_prep (&L, &U, A, msg)) ;
            GRB_TRY (GrB_mxm (C, L, NULL, semiring, L, U, GrB_DESC_ST1)) ;
            GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case LAGraph_TriangleCount_SandiaDot2: // 6: sum (sum ((U * L') .* U))

            // using the masked dot product
            LG_TRY (tricount_prep (&L, &U, A, msg)) ;
            GRB_TRY (GrB_mxm (C, U, NULL, semiring, U, L, GrB_DESC_ST1)) ;
            GRB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    (*ntriangles) = (uint64_t) ntri ;
    return (GrB_SUCCESS) ;
}
