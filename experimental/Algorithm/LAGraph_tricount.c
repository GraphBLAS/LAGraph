//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2020 LAGraph Contributors.

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

// LAGraph_tricount: count the number of triangles in a graph,
// Contributed by Tim Davis, Texas A&M.

// Given a symmetric graph A with no-self edges, LAGraph_tricount counts the
// number of triangles in the graph.  A triangle is a clique of size three,
// that is, 3 nodes that are all pairwise connected.

// One of 6 methods are used, defined below where L and U are the strictly
// lower and strictly upper triangular parts of the symmetrix matrix A,
// respectively.  Each method computes the same result, ntri:

//  1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia:     ntri = sum (sum ((L * L) .* L))
//  4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
//  5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  SandiaDot2: ntri = sum (sum ((U * L') .* U)).  Note that U=L'.

// A is a square symmetric matrix, of any type.  Its values are ignored,
// (assuming v3.2.0 of SuiteSparse:GraphBLAS is used); otherwise, A must be
// binary.  Results are undefined for methods 1 and 2 if self-edges exist in A.
// Results are undefined for all methods if A is unsymmetric.

// TODO use an enum for the above methods.

// All matrices are assumed to be in CSR format (GxB_BY_ROW in
// SuiteSparse:GraphBLAS).  The 6 methods work fine if the matrices are in CSC
// format; just the underlying algorithms employed inside SuiteSparse:GraphBLAS
// will differ (dot product vs saxpy, for example).  If L and U are in CSC
// format, then the "Dot" methods would use an outer product approach, which is
// slow in SuiteSparse:GraphBLAS (requiring an explicit transpose).  The
// auto-sort rule probably needs to be reversed, if A is in CSC format (this is
// not yet tested).

// Methods 1 and 2 are much slower than methods 3 to 6 and take more memory.
// Methods 3 to 6 take a little less memory than methods 1 and 2, are by far
// the fastest methods in general.  The methods 3 and 5 compute the same
// intermediate matrix (L*L), and differ only in the way the matrix
// multiplication is done.  Method 3 uses an outer-product method (Gustavson's
// method).  Method 5 uses dot products (assuming both matrices are in CSR
// format) and does not explicitly transpose U.  They are called the "Sandia"
// method since matrices in the KokkosKernels are stored in compressed-sparse
// row form, so (L*L).*L in the KokkosKernel method is equivalent to (L*L).*L
// in SuiteSparse:GraphBLAS when the matrices in SuiteSparse:GraphBLAS are in
// their default format (also by row).

// The new GxB_PAIR_INT64 binary operator in SuiteSparse:GraphBLAS v3.2.0 is
// used in the semiring, if available.  This is the function f(x,y)=1, so the
// values of A are not accessed.  They can have any values and any type.  Only
// the structure of A.  Otherwise, without this operator, the input matrix A
// must be binary.

// Reference: Wolf, Deveci, Berry, Hammond, Rajamanickam, 'Fast linear algebra-
// based triangle counting with KokkosKernels', IEEE HPEC'17,
// https://dx.doi.org/10.1109/HPEC.2017.8091043,

#include "LAGraph_internal.h"
#include "GB_msort_2.h"

//------------------------------------------------------------------------------
// tricount_prep: construct L and U
//------------------------------------------------------------------------------

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL        \
    GrB_free (&thunk) ;         \
    GrB_free (L) ;              \
    GrB_free (U) ;

static GrB_Info tricount_prep
(
    GrB_Matrix *L,
    GrB_Matrix *U,
    GrB_Matrix A
)
{
    GrB_Index n, *I = NULL, *J = NULL ;
    bool *X = NULL ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,0,1) )

        //----------------------------------------------------------------------
        // build L and/or U with GxB_select
        //----------------------------------------------------------------------

        GxB_Scalar thunk ;
        LAGr_Matrix_nrows (&n, A) ;
        LAGr_Scalar_new (&thunk, GrB_INT64) ;

        if (L != NULL)
        {
            // L = tril (A,-1)
            LAGr_Matrix_new (L, GrB_BOOL, n, n) ;
            LAGr_Scalar_setElement (thunk, -1) ;
            LAGr_select (*L, NULL, NULL, GxB_TRIL, A, thunk, NULL) ;
        }

        if (U != NULL)
        {
            // U = triu (A,1)
            LAGr_Matrix_new (U, GrB_BOOL, n, n) ;
            LAGr_Scalar_setElement (thunk, 1) ;
            LAGr_select (*U, NULL, NULL, GxB_TRIU, A, thunk, NULL) ;
        }

        LAGr_free (&thunk) ;

    #else

        //----------------------------------------------------------------------
        // build L and U with extractTuples (slower than GxB_select)
        //----------------------------------------------------------------------

        GrB_Vector thunk ;
        LAGr_Matrix_nrows (&n, A) ;
        if (L != NULL || U != NULL)
        {
            GrB_Index nvals ;
            LAGr_Matrix_nvals (&nvals, A) ;
            I = LAGraph_malloc (nvals, sizeof (GrB_Index)) ;
            J = LAGraph_malloc (nvals, sizeof (GrB_Index)) ;
            X = LAGraph_malloc (nvals, sizeof (bool)) ;
            if (I == NULL || J == NULL || X == NULL)
            {
                LAGRAPH_ERROR ("out of memory") ;
            }

            LAGr_Matrix_extractTuples (I, J, X, &nvals, A) ;

            // remove entries in the upper triangular part
            nedges = 0 ;
            for (int64_t k = 0 ; k < nvals ; k++)
            {
                if (I [k] > J [k])
                {
                    // keep this entry
                    I [nedges] = I [k] ;
                    J [nedges] = J [k] ;
                    X [nedges] = X [k] ;
                    nedges++ ;
                }
            }

            if (L != NULL)
            {
                LAGr_Matrix_new (L, GrB_BOOL, n, n) ;
                LAGr_Matrix_build (*L, I, J, X, nedges, GrB_LOR) ;
            }

            if (U != NULL)
            {
                LAGr_Matrix_new (U, GrB_BOOL, n, n) ;
                LAGr_Matrix_build (*U, J, I, X, nedges, GrB_LOR) ;
            }

            LAGRAPH_FREE (I) ;
            LAGRAPH_FREE (J) ;
            LAGRAPH_FREE (X) ;
        }
    #endif
}

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL        \
    GrB_free (&C) ;             \
    GrB_free (&L) ;             \
    GrB_free (&T) ;             \
    GrB_free (&U) ;             \
    LAGRAPH_FREE (W0) ;         \
    LAGRAPH_FREE (W1) ;         \
    LAGRAPH_FREE (P) ;          \
    LAGRAPH_FREE (D) ; 

#if 0
// easy mode:
LAGr_info LAGraph_tricount 
(
    uint64_t *ntriangles,   // # of triangles
    LAGr_Graph G,           // a graph
    LAGr_descriptor d
) ;

LAGr_info LAGraph_tricount 
(
    uint64_t *ntriangles,   // # of triangles
    bool directed,
    LAGr_Matrix A           // adj matrix of an directed graph
) ;
#endif

GrB_Info LAGraph_tricount   // count # of triangles
(
    int64_t *ntri,          // # of triangles
    const int method,       // 1 to 6, see above
    int sorting,            //  0: no sort
                            //  1: sort by degree, ascending order
                            // -1: sort by degree, descending order
                            //  2: auto selection: no sort if rule is not
                            // triggered.  Otherise: sort in ascending order
                            // for methods 3 and 5, descending ordering for
                            // methods 4 and 6.
    const int64_t *degree,  // degree of each node, may be NULL if sorting==0.
                            // of size n, unmodified. 
    const GrB_Matrix A_in   // input matrix, must be symmetric, no diag entries
)
{

    //--------------------------------------------------------------------------
    // check inputs and initialize
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Index n ;
    GrB_Matrix C = NULL, L = NULL, U = NULL, T = NULL, A = NULL ;
    int64_t *P = NULL, *D = NULL, *W0 = NULL, *W1 = NULL ;
    LAGr_Matrix_nrows (&n, A_in) ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,2,0) )
        // the PAIR function is f(x,y)=1, ignoring input values and type
        GrB_Semiring   semiring = GxB_PLUS_PAIR_INT64 ;
    #else
        // f(x,y)=x*y, so x and y must be 1 to compute the correct count, and
        // thus the input matrix A must be binary.
        GrB_Semiring   semiring = LAGraph_PLUS_TIMES_INT64 ;
    #endif

    GrB_Monoid sum = LAGraph_PLUS_INT64_MONOID ;
    LAGr_Matrix_new (&C, GrB_INT64, n, n) ;

    //--------------------------------------------------------------------------
    // heuristic sort rule
    //--------------------------------------------------------------------------

    if (sorting == 2)
    {
        // auto selection of sorting method
        sorting = 0 ;       // default is not to sort

        if (method >= 3 && method <= 6)
        {
            // This rule is very similar to Scott Beamer's rule in the GAP TC
            // benchmark, except that it is extended to handle the ascending
            // sort needed by methods 3 and 5.  It also uses a stricter rule,
            // since the performance of triangle counting in GraphBLAS is less
            // sensitive to the sorting as compared to the GAP algorithm.  This
            // is because the dot products in GraphBLAS use binary search if
            // one vector is very sparse compared to the other.  As a result,
            // GraphBLAS needs the sort for fewer matrices, as compared to the
            // GAP algorithm.

            // With this rule, the GAP-kron and GAP-twitter matrices are
            // sorted, and the others remain unsorted.  With the rule in the
            // GAP tc.cc benchmark, GAP-web is also sorted, but it is not
            // sorted here.

            #define NSAMPLES 1000
            GrB_Index nvals ;
            LAGr_Matrix_nvals (&nvals, A_in) ;
            if (n > NSAMPLES && ((double) nvals / ((double) n)) >= 10)
            {
                // pick 1000 nodes at random and determine their degree
                // struct drand48_data buffer ;
                // srand48_r ((long int) n, &buffer) ;
                uint64_t seed = n ;
                int64_t samples [NSAMPLES] ;
                int64_t dsum = 0 ;
                for (int k = 0 ; k < NSAMPLES ; k++)
                {
                    uint64_t result = LAGraph_rand64 (&seed) ;
                    // lrand48_r (&buffer, &result) ;
                    int64_t i = result % n ;
                    int64_t d = degree [i] ;
                    samples [k] = d ;
                    dsum += d ;
                }

                // find the average degree
                double sample_average = ((double) dsum) / NSAMPLES ;

                // find the median degree
                GB_qsort_1a (samples, NSAMPLES) ;
                double sample_median = (double) samples [NSAMPLES/2] ;

                printf ("average degree: %g\n", sample_average) ;
                printf ("median  degree: %g\n", sample_median) ;

                // sort if the average degree is very high compared to the
                // median
                if (sample_average > 4 * sample_median)
                {
                    switch (method)
                    {
                        case 3:  sorting =  1 ; break ;     // sort ascending
                        case 4:  sorting = -1 ; break ;     // sort descending
                        case 5:  sorting =  1 ; break ;     // sort ascending
                        case 6:  sorting = -1 ; break ;     // sort descending
                        default: sorting =  0 ; break ;     // no sort
                    }
                }
            }
        }

        printf ("auto sorting: %d: ", sorting) ;
        if (sorting == 0) printf ("none") ;
        else if (sorting == -1) printf ("descending") ;
        else if (sorting ==  1) printf ("ascending") ;
        printf ("\n") ;
    }

    //--------------------------------------------------------------------------
    // sort the input matrix, if requested
    //--------------------------------------------------------------------------

    if (sorting != 0)
    {

        // decide how many threads to use
        #define CHUNK (64*1024)
        int nthreads = LAGraph_get_nthreads ( ) ;
        nthreads = LAGRAPH_MIN (nthreads, n/CHUNK) ;
        nthreads = LAGRAPH_MAX (nthreads, 1) ;

        // allocate workspace
        P  = LAGraph_malloc (n, sizeof (int64_t)) ;
        D  = LAGraph_malloc (n, sizeof (int64_t)) ;
        W0 = LAGraph_malloc (n, sizeof (int64_t)) ;
        W1 = LAGraph_malloc (n, sizeof (int64_t)) ;
        if (P == NULL || D == NULL || W0 == NULL || W1 == NULL)
        {
            // out of memory
            LAGRAPH_FREE_ALL ;
            return (GrB_OUT_OF_MEMORY) ;
        }

        // construct the pair [D,P] to sort
        if (sorting > 0)
        {
            printf ("sort ascending\n") ;
            // sort [D,P] in ascending order of degree, tie-breaking on P
            #pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int64_t k = 0 ; k < n ; k++)
            {
                D [k] = degree [k] ;
                P [k] = k ;
            }
        }
        else
        {
            printf ("sort descending\n") ;
            // sort [D,P] in descending order of degree, tie-breaking on P
            #pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int64_t k = 0 ; k < n ; k++)
            {
                D [k] = -degree [k] ;
                P [k] = k ;
            }
        }

        GB_msort_2 (D, P, W0, W1, n, nthreads) ;

        // T = A_in (P,P) and typecast to boolean
        LAGr_Matrix_new (&T, GrB_BOOL, n, n) ;
        LAGr_extract (T, NULL, NULL, A_in, P, n, P, n, NULL) ;
        A = T ;
    }
    else
    {
        // use the input matrix as-is
        A = A_in ;
    }

    // free workspace
    LAGRAPH_FREE (W0) ; 
    LAGRAPH_FREE (W1) ; 
    LAGRAPH_FREE (D) ; 
    LAGRAPH_FREE (P) ; 

    //--------------------------------------------------------------------------
    // count triangles
    //--------------------------------------------------------------------------

    switch (method)
    {
        #if 0
        // case 0:  // minitri:    ntri = nnz (A*E == 2) / 3

            // This method requires the incidence matrix E.  It is very slow
            // compared to the other methods.  The construction of E was done
            // in the Test/Tricount/*.c driver, and it hasn't been added here.
            LAGr_Matrix_ncols (&ne, E) ;
            LAGr_free (&C) ;
            LAGr_Matrix_new (&C, GrB_INT64, n, ne) ;
            LAGr_mxm (C, NULL, NULL, semiring, A, E, NULL) ;
            LAGr_Matrix_new (&S, GrB_BOOL, n, ne) ;
            LAGr_apply (S, NULL, NULL, LAGraph_ISTWO_INT64, C, NULL) ;
            LAGr_reduce (ntri, NULL, sum, S, NULL) ;
            (*ntri) /= 3 ;
            break ;
        #endif

        case 1:  // Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6

            LAGr_mxm (C, A, NULL, semiring, A, A, GrB_DESC_S) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            (*ntri) /= 6 ;
            break ;

        case 2:  // Cohen:      ntri = sum (sum ((L * U) .* A)) / 2

            LAGRAPH_OK (tricount_prep (&L, &U, A)) ;
            LAGr_mxm (C, A, NULL, semiring, L, U, GrB_DESC_S) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            (*ntri) /= 2 ;
            break ;

        case 3:  // Sandia:     ntri = sum (sum ((L * L) .* L))

            // using the masked saxpy3 method
            LAGRAPH_OK (tricount_prep (&L, NULL, A)) ;
            LAGr_mxm (C, L, NULL, semiring, L, L, GrB_DESC_S) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        case 4:  // Sandia2:    ntri = sum (sum ((U * U) .* U))

            // using the masked saxpy3 method
            LAGRAPH_OK (tricount_prep (NULL, &U, A)) ;
            LAGr_mxm (C, U, NULL, semiring, U, U, GrB_DESC_S) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        case 5:  // SandiaDot:  ntri = sum (sum ((L * U') .* L))

            // This tends to be the fastest method, for most matrices, but the
            // Dot2 method is also very fast.

            // using the masked dot product
            LAGRAPH_OK (tricount_prep (&L, &U, A)) ;
            LAGr_mxm (C, L, NULL, semiring, L, U, GrB_DESC_ST1) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        case 6:  // SandiaDot2: ntri = sum (sum ((U * L') .* U))

            // using the masked dot product
            LAGRAPH_OK (tricount_prep (&L, &U, A)) ;
            LAGr_mxm (C, U, NULL, semiring, U, L, GrB_DESC_ST1) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        default:    // invalid method

            LAGRAPH_FREE_ALL ;
            return (GrB_INVALID_VALUE) ;
            break ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

