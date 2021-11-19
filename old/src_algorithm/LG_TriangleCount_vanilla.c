//------------------------------------------------------------------------------
// LG_TriangleCount_vanilla: count the number of triangles in a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_tricount: count the number of triangles in a graph,
// Contributed by Scott McMillan, CMU and Tim Davis, Texas A&M.

// Given a symmetric graph A with no-self edges, LAGraph_TriangleCount_methods
// counts the number of triangles in the graph.  A triangle is a clique of size
// three, that is, 3 nodes that are all pairwise connected.

// One of 6 methods are used, defined below where L and U are the strictly
// lower and strictly upper triangular parts of the symmetrix matrix A,
// respectively.  Each method computes the same result, ntri:

//  0:  minitri:    ntri = nnz (A*E == 2) / 3 ; this method is disabled.
//  1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia:     ntri = sum (sum ((L * L) .* L))
//  4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
//  5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  SandiaDot2: ntri = sum (sum ((U * L') .* U)).  Note that U=L'.

// A is a square symmetric matrix, of any type.  Its values are ignored.
// Results are undefined for methods 1 and 2 if self-edges exist in A.  Results
// are undefined for all methods if A is unsymmetric.  Method 0 (minitri) is
// not yet available, since it requires G to include both an adjacency matrix
// and an incidence matrix (in any case, minitri is the slowest method and
// is included only for reference).

// The Sandia* methods all tend to be faster than the Burkhard or Cohen
// methods.  For the largest graphs, SandiaDot tends to be fastest, except for
// the GAP-urand matrix, where the saxpy-based Sandia method (L*L.*L) is
// fastest.  For many small graphs, the saxpy-based Sandia and Sandia2 methods
// are often faster that the dot-product-base methods.

// TODO use an enum for the above methods.

// Reference (for the "Sandia*" methods): Wolf, Deveci, Berry, Hammond,
// Rajamanickam, 'Fast linear algebra- based triangle counting with
// KokkosKernels', IEEE HPEC'17, https://dx.doi.org/10.1109/HPEC.2017.8091043,

#define LAGraph_FREE_ALL        \
{                               \
    GrB_free (L) ;              \
    GrB_free (U) ;              \
}

#include <LG_internal.h>

// FIXME set NO_GRB_SELECT according to GRB_VERSION
#define NO_GRB_SELECT 1

//------------------------------------------------------------------------------
// tricount_prep: construct L and U for LAGraph_TriangleCount_Methods
//------------------------------------------------------------------------------

static int tricount_prep_vanilla        // return 0 if successful, < 0 on error
(
    GrB_Matrix *L,      // if present, compute L = tril (A,-1)
    GrB_Matrix *U,      // if present, compute U = triu (A, 1)
    GrB_Matrix  A,      // input matrix
    GrB_Type    type,
    char *msg
)
{
    if ((L == NULL) && (U == NULL))
    {
        return -1;
    }

    GrB_Index n, nvals;
    GrB_TRY( GrB_Matrix_nrows (&n, A) );
    GrB_TRY( GrB_Matrix_nvals (&nvals, A) );

#if defined(NO_GRB_SELECT)
    // Allocate space to do the extract
    size_t     sz  = 0;
    size_t     vsz = 0;
    GrB_Index *row_ids;
    GrB_Index *col_ids;
    GrB_Index  nv   = nvals;
    row_ids = (GrB_Index*) malloc(nvals*sizeof(GrB_Index));
    col_ids = (GrB_Index*) malloc(nvals*sizeof(GrB_Index));

    if (type == GrB_BOOL)
    {
        bool *values = (bool*)malloc(nvals*sizeof(bool));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_INT8)
    {
        int8_t *values = malloc(nvals*sizeof(int8_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_INT16)
    {
        int16_t *values = malloc(nvals*sizeof(int16_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_INT32)
    {
        int32_t *values = malloc(nvals*sizeof(int32_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_INT64)
    {
        int64_t *values = malloc(nvals*sizeof(int64_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_UINT8)
    {
        uint8_t *values = malloc(nvals*sizeof(uint8_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_UINT16)
    {
        uint16_t *values = malloc(nvals*sizeof(uint16_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_UINT32)
    {
        uint32_t *values = malloc(nvals*sizeof(uint32_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_UINT64)
    {
        uint64_t *values = malloc(nvals*sizeof(uint64_t));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_FP32)
    {
        float *values = malloc(nvals*sizeof(float));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else if (type == GrB_FP64)
    {
        double *values = malloc(nvals*sizeof(double));
        GrB_TRY( GrB_Matrix_extractTuples(row_ids, col_ids, values, &nv, A) );
        free(values);
    }
    else
    {
        // unknown type
        free(row_ids);
        free(col_ids);
        return -2;
    }

    GrB_Index l_nvals = 0;
    GrB_Index u_nvals = 0;
    for (GrB_Index ix = 0; ix < nvals; ++ix)
    {
        if (row_ids[ix] > col_ids[ix])
        {
            ++l_nvals;
        }
        else if (row_ids[ix] < col_ids[ix])
        {
            ++u_nvals;
        }
    }
#endif

    if (L != NULL)
    {
        // L = tril (A,-1)
        GrB_TRY( GrB_Matrix_new (L, GrB_BOOL, n, n) );
#if defined(NO_GRB_SELECT)
        size_t     lsz  = 0;
        size_t     lvsz = 0;
        GrB_Index *lrow_ids;
        GrB_Index *lcol_ids;
        bool      *lvals;
        lrow_ids = (GrB_Index*) malloc(l_nvals*sizeof(GrB_Index));
        lcol_ids = (GrB_Index*) malloc(l_nvals*sizeof(GrB_Index));
        lvals    = (bool*)      malloc(l_nvals*sizeof(bool));
        GrB_Index idx = 0;
        for (GrB_Index ix = 0; ix < nvals; ++ix)
        {
            if (row_ids[ix] > col_ids[ix])
            {
                lrow_ids[idx] = row_ids[ix];
                lcol_ids[idx] = col_ids[ix];
                lvals[idx] = true;
                ++idx;
            }
        }
        GrB_TRY( GrB_Matrix_build(*L, lrow_ids, lcol_ids, lvals, idx,
                                  GrB_SECOND_BOOL) );
        free(lrow_ids);
        free(lcol_ids);
        free(lvals);
#else
        GrB_TRY( GrB_select (*L, GrB_NULL, GrB_NULL,
                             GrB_TRIL_BOOL, A, 0UL, GrB_NULL) );
#endif
    }

    if (U != NULL)
    {
        // U = triu (A,1)
        GrB_TRY( GrB_Matrix_new (U, GrB_BOOL, n, n) );
#if defined(NO_GRB_SELECT)
        size_t     usz  = 0;
        size_t     uvsz = 0;
        GrB_Index *urow_ids;
        GrB_Index *ucol_ids;
        bool      *uvals;
        urow_ids = (GrB_Index*) malloc(u_nvals*sizeof(GrB_Index));
        ucol_ids = (GrB_Index*) malloc(u_nvals*sizeof(GrB_Index));
        uvals    = (bool*)      malloc(u_nvals*sizeof(bool));
        GrB_Index idx = 0;
        for (GrB_Index ix = 0; ix < nvals; ++ix)
        {
            if (row_ids[ix] < col_ids[ix])
            {
                urow_ids[idx] = row_ids[ix];
                ucol_ids[idx] = col_ids[ix];
                uvals[idx] = true;
                ++idx;
            }
        }
        GrB_TRY( GrB_Matrix_build(*U, urow_ids, ucol_ids, uvals, idx,
                                  GrB_SECOND_BOOL) );
        free(urow_ids);
        free(ucol_ids);
        free(uvals);
#else
        GrB_TRY( GrB_select (*U, GrB_NULL, GrB_NULL,
                             GrB_TRIU_BOOL, A, 0UL, GrB_NULL) );
#endif
    }

#if defined(NO_GRB_SELECT)
    free(row_ids);
    free(col_ids);
#endif

    return (0) ;
}

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL                    \
{                                           \
    GrB_free (&C) ;                         \
    GrB_free (&L) ;                         \
    GrB_free (&T) ;                         \
    GrB_free (&U) ;                         \
    LAGraph_Free ((void **) &P) ;   \
}

int LG_TriangleCount_vanilla   // returns 0 if successful, < 0 if failure
(
    uint64_t *ntriangles,   // # of triangles
    // input:
    LAGraph_Graph G,
    int method,             // selects the method to use (TODO: enum)
    // input/output:
    int *presort,           // controls the presort of the graph (TODO: enum)
        //  0: no sort
        //  1: sort by degree, ascending order
        // -1: sort by degree, descending order
        //  2: auto selection: no sort if rule is not triggered.  Otherise:
        //  sort in ascending order for methods 3 and 5, descending ordering
        //  for methods 4 and 6.  On output, presort is modified to reflect the
        //  sorting method used (0, -1, or 1).  If presort is NULL on input, no
        //  sort is performed.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL, A ;
    int64_t *P = NULL ;
    LG_CHECK ((method < 1) || (method > 6), -101, "method is invalid");
    LG_CHECK (LAGraph_CheckGraph (G, msg), -102, "graph is invalid") ;
    LG_CHECK (ntriangles == NULL, -103, "ntriangles is null") ;
    LG_CHECK (G->ndiag != 0, -104, "G->ndiag must be zero") ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_structure_is_symmetric == LAGRAPH_TRUE))
    {
        // the structure of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -105, "G->A must be symmetric") ;
    }

    GrB_Vector Degree = G->rowdegree ;
    bool auto_sort = (presort != NULL && (*presort) == 2) ;
    if (auto_sort && method >= 3 && method <= 6)
    {
        LG_CHECK (Degree == NULL, -106, "G->rowdegree must be defined") ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_new (&C, GrB_INT64, n, n)) ;
    // FIXME create the PLUS_ONEB_INT64 for v2.0 C API
    GrB_Semiring semiring = GrB_PLUS_TIMES_SEMIRING_INT64 ;
    GrB_Monoid monoid = GrB_PLUS_MONOID_INT64 ;

    //--------------------------------------------------------------------------
    // heuristic sort rule
    //--------------------------------------------------------------------------

    if (auto_sort)
    {
        // auto selection of sorting method
        (*presort) = 0 ;       // default is not to sort

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
            // GAP tc.cc benchmark, GAP-kron and GAP-twitter are sorted, and so
            // is GAP-web, but GAP-web is not sorted here.

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
                        case 3:  (*presort) =  1 ; break ;  // sort ascending
                        case 4:  (*presort) = -1 ; break ;  // sort descending
                        case 5:  (*presort) =  1 ; break ;  // sort ascending
                        case 6:  (*presort) = -1 ; break ;  // sort descending
                        default: (*presort) =  0 ; break ;  // no sort
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // sort the input matrix, if requested
    //--------------------------------------------------------------------------

    if (presort != NULL && (*presort) != 0)
    {
        // P = permutation that sorts the rows by their degree
        LAGraph_TRY (LAGraph_SortByDegree (&P, G, true, (*presort) > 0, msg)) ;

        // T = A (P,P) and typecast to boolean
        GrB_TRY (GrB_Matrix_new (&T, GrB_BOOL, n, n)) ;
        GrB_TRY (GrB_extract (T, NULL, NULL, A, (GrB_Index *) P, n,
            (GrB_Index *) P, n, NULL)) ;
        A = T ;

        // free workspace
        LAGraph_Free((void **) &P) ;
        P = NULL;
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

            LAGraph_TRY (tricount_prep_vanilla (&L, &U, A, G->A_type, msg)) ;
            GrB_TRY (GrB_mxm (C, A, NULL, semiring, L, U, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            ntri /= 2 ;
            break ;

        case 3:  // Sandia:     ntri = sum (sum ((L * L) .* L))

            // using the masked saxpy3 method
            LAGraph_TRY (tricount_prep_vanilla (&L, NULL, A, G->A_type, msg)) ;
            GrB_TRY (GrB_mxm (C, L, NULL, semiring, L, L, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case 4:  // Sandia2:    ntri = sum (sum ((U * U) .* U))

            // using the masked saxpy3 method
            LAGraph_TRY (tricount_prep_vanilla (NULL, &U, A, G->A_type, msg)) ;
            GrB_TRY (GrB_mxm (C, U, NULL, semiring, U, U, GrB_DESC_S)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case 5:  // SandiaDot:  ntri = sum (sum ((L * U') .* L))

            // This tends to be the fastest method for most large matrices, but
            // the SandiaDot2 method is also very fast.

            // using the masked dot product
            LAGraph_TRY (tricount_prep_vanilla (&L, &U, A, G->A_type, msg)) ;
            GrB_TRY (GrB_mxm (C, L, NULL, semiring, L, U, GrB_DESC_ST1)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        case 6:  // SandiaDot2: ntri = sum (sum ((U * L') .* U))

            // using the masked dot product
            LAGraph_TRY (tricount_prep_vanilla (&L, &U, A, G->A_type, msg)) ;
            GrB_TRY (GrB_mxm (C, U, NULL, semiring, U, L, GrB_DESC_ST1)) ;
            GrB_TRY (GrB_reduce (&ntri, NULL, monoid, C, NULL)) ;
            break ;

        default:    // invalid method

            LG_CHECK (false, -101, "invalid method") ;
            break ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_ALL ;
    (*ntriangles) = (uint64_t) ntri ;
    return (0) ;
}
