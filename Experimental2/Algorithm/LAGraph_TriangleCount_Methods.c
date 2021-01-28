//------------------------------------------------------------------------------
// LAGraph_TriangleCount_Methods: count the number of triangles in a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// LAGraph_tricount: count the number of triangles in a graph,
// Contributed by Tim Davis, Texas A&M.

// Given a symmetric graph A with no-self edges, LAGraph_TriangleCount_methods
// counts the number of triangles in the graph.  A triangle is a clique of size
// three, that is, 3 nodes that are all pairwise connected.

// One of TODO methods are used, defined below where L and U are the strictly
// lower and strictly upper triangular parts of the symmetrix matrix A,
// respectively.  Each method computes the same result, ntri:

//  1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia:     ntri = sum (sum ((L * L) .* L))
//  4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
//  5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  SandiaDot2: ntri = sum (sum ((U * L') .* U)).  Note that U=L'.

// A is a square symmetric matrix, of any type.  Its values are ignored.
// Results are undefined for methods 1 and 2 if self-edges exist in A.
// Results are undefined for all methods if A is unsymmetric.

// TODO use an enum for the above methods.

// Reference (for the "Sandia*" methods): Wolf, Deveci, Berry, Hammond,
// Rajamanickam, 'Fast linear algebra- based triangle counting with
// KokkosKernels', IEEE HPEC'17, https://dx.doi.org/10.1109/HPEC.2017.8091043,

#include "LG_internal.h"

//------------------------------------------------------------------------------
// tricount_prep: construct L and U for LAGraph_TriangleCount_Methods
//------------------------------------------------------------------------------

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL        \
{                               \
    GrB_free (&thunk) ;         \
    GrB_free (L) ;              \
    GrB_free (U) ;              \
}

static int tricount_prep
(
    GrB_Matrix *L,
    GrB_Matrix *U,
    GrB_Matrix A,
    char *msg
)
{
    GrB_Index n ;
    GxB_Scalar thunk ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GxB_Scalar_new (&thunk, GrB_INT64)) ;

    if (L != NULL)
    {
        // L = tril (A,-1)
        GrB_TRY (GrB_Matrix_new (L, GrB_BOOL, n, n)) ;
        GrB_TRY (GxB_Scalar_setElement (thunk, -1)) ;
        GrB_TRY (GxB_select (*L, NULL, NULL, GxB_TRIL, A, thunk, NULL)) ;
    }

    if (U != NULL)
    {
        // U = triu (A,1)
        GrB_TRY (GrB_Matrix_new (U, GrB_BOOL, n, n)) ;
        GrB_TRY (GxB_Scalar_setElement (thunk, 1)) ;
        GrB_TRY (GxB_select (*U, NULL, NULL, GxB_TRIU, A, thunk, NULL)) ;
    }

    GrB_free (&thunk) ;
    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL        \
{                               \
    GrB_free (&C) ;             \
    GrB_free (&L) ;             \
    GrB_free (&T) ;             \
    GrB_free (&U) ;             \
}

int LAGraph_TriangleCount_Methods   // returns -1 on failure, 0 if successful
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph G,
    int method,             // selects the method to use (TODO: enum)
    int presort,            // controls the presort of the graph (TODO: enum)
        //  0: no sort
        //  1: sort by degree, ascending order
        // -1: sort by degree, descending order
        //  2: auto selection: no sort if rule is not triggered.  Otherise:
        //  sort in ascending order for methods 3 and 5, descending ordering
        //  for methods 4 and 6.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL, A ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;

    LAGraph_DisplayGraph (G, 2, msg) ;
    printf ("kind %d\n", G->kind) ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_pattern_is_symmetric == LAGRAPH_TRUE))
    {
        // the pattern of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -1, "adjacency matrix must be symmetric") ;
    }

    GrB_Vector Degree = G->rowdegree ;
    if (presort == 2 && method >= 3 && method <= 6)
    {
        LG_CHECK (Degree == NULL, -1, "G->rowdegree must be defined") ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_Semiring semiring = GxB_PLUS_PAIR_INT64 ;
    GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;
    GrB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;

    //--------------------------------------------------------------------------
    // heuristic sort rule
    //--------------------------------------------------------------------------

    if (presort == 2)
    {
        // auto selection of sorting method
        presort = 0 ;       // default is not to sort

        if (method >= 3 && method <= 6)
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
            // GAP tc.cc benchmark, GAP-web is also sorted, but it is not
            // sorted here.

            #define NSAMPLES 1000
            GrB_Index nvals ;
            GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
            if (n > NSAMPLES && ((double) nvals / ((double) n)) >= 10)
            {
                // estimate the mean and median degrees
                double mean, median ;
                LAGraph_TRY (LAGraph_SampleDegree (&mean, &median,
                    G, true, NSAMPLES, n, msg)) ;
                // sort if the average degree is very high vs the median
                if (mean > 4 * median)
                {
                    switch (method)
                    {
                        case 3:  presort =  1 ; break ;     // sort ascending
                        case 4:  presort = -1 ; break ;     // sort descending
                        case 5:  presort =  1 ; break ;     // sort ascending
                        case 6:  presort = -1 ; break ;     // sort descending
                        default: presort =  0 ; break ;     // no sort
                    }
                }
            }
        }

        // printf ("auto sorting: %d: ", presort) ;
        // if (presort == 0) printf ("none") ;
        // else if (presort == -1) printf ("descending") ;
        // else if (presort ==  1) printf ("ascending") ;
        // printf ("\n") ;
    }

    //--------------------------------------------------------------------------
    // sort the input matrix, if requested
    //--------------------------------------------------------------------------

    if (presort != 0)
    {
        // P = permutation that sorts the rows by their degree
        int64_t *P ;
        LAGraph_TRY (LAGraph_SortByDegree (&P, G, true, presort > 0, msg)) ;

        // T = A (P,P) and typecast to boolean
        GrB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;
        GrB_TRY (GrB_extract (T, NULL, NULL, A, P, n, P, n, NULL)) ;
        A = T ;

        // free workspace
        LAGraph_FREE (P) ; 
    }

    //--------------------------------------------------------------------------
    // count triangles
    //--------------------------------------------------------------------------

    int64_t ntri ;

    switch (method)
    {
        #if 0
        // case 0:  // minitri:    ntri = nnz (A*E == 2) / 3

            // This method requires the incidence matrix E.  It is very slow
            // compared to the other methods.  The LAGraph_Graph does not yet
            // include an incidence matrix, so this method is here only for
            // reference and possible future use.
            GrB_TRY (GrB_Matrix_ncols (&ne, E)) ;
            GrB_TRY (GrB_free (&C)) ;
            GrB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, ne)) ;
            GrB_TRY (GrB_mxm (C, NULL, NULL, semiring, A, E, NULL)) ;
            GrB_TRY (GrB_Matrix_new (&S, GrB_BOOL, n, ne)) ;
            GrB_TRY (GrB_apply (S, NULL, NULL, LAGraph_ISTWO_INT64, C, NULL)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, S, NULL)) ;
            ntri /= 3 ;
            break ;
        #endif

        case 1:  // Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6

            GrB_TRY (GrB_mxm (C, A, NULL, semiring, A, A, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            ntri /= 6 ;
            break ;

        case 2:  // Cohen:      ntri = sum (sum ((L * U) .* A)) / 2

            LAGraph_TRY (tricount_prep (&L, &U, A, msg)) ;
            GrB_TRY (GrB_mxm (C, A, NULL, semiring, L, U, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            ntri /= 2 ;
            break ;

        case 3:  // Sandia:     ntri = sum (sum ((L * L) .* L))

            // using the masked saxpy3 method
            LAGraph_TRY (tricount_prep (&L, NULL, A, msg)) ;
            GrB_TRY (GrB_mxm (C, L, NULL, semiring, L, L, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case 4:  // Sandia2:    ntri = sum (sum ((U * U) .* U))

            // using the masked saxpy3 method
            LAGraph_TRY (tricount_prep (NULL, &U, A, msg)) ;
            GrB_TRY (GrB_mxm (C, U, NULL, semiring, U, U, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case 5:  // SandiaDot:  ntri = sum (sum ((L * U') .* L))

            // This tends to be the fastest method for most large matrices, but
            // the SandiaDot2 method is also very fast.

            // using the masked dot product
            LAGraph_TRY (tricount_prep (&L, &U, A, msg)) ;
            GrB_TRY (GrB_mxm (C, L, NULL, semiring, L, U, GrB_DESC_ST1)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case 6:  // SandiaDot2: ntri = sum (sum ((U * L') .* U))

            // using the masked dot product
            LAGraph_TRY (tricount_prep (&L, &U, A, msg)) ;
            GrB_TRY (GrB_mxm (C, U, NULL, semiring, U, L, GrB_DESC_ST1)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        default:    // invalid method

            LG_CHECK (false, -1, "invalid method") ;
            break ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    (*ntriangles) = (uint64_t) ntri ;
    return (0) ;
}

