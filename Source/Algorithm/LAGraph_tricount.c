//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
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

// LAGraph_tricount: count the number of triangles in a graph,
// Contributed by Tim Davis, Texas A&M.

// Given a symmetric binary graph A with no-self edges, LAGraph_tricount counts
// the exact number of triangles in the graph.  A triangle is a clique of size
// three, that is, 3 nodes that are all pairwise connected.

// On input, the L and U matrices are the strictly lower and strictly upper
// triangular parts of the symmetrix matrix A, respectively.

// One of 6 methods are used.  Each computes the same result, ntri:

//  0:  minitri:    ntri = nnz (A*E == 2) / 3
//  1:  Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6
//  2:  Cohen:      ntri = sum (sum ((L * U) .* A)) / 2
//  3:  Sandia:     ntri = sum (sum ((L * L) .* L))
//  4:  Sandia2:    ntri = sum (sum ((U * U) .* U))
//  5:  SandiaDot:  ntri = sum (sum ((L * U') .* L)).  Note that L=U'.
//  6:  SandiaDot2: ntri = sum (sum ((U * L') .* U))

// TODO use an enum for the above methods.

// All matrices are assumed to be in CSR format (GxB_BY_ROW in
// SuiteSparse:GraphBLAS).  The 6 methods work fine if the matrices are in CSC
// format; just the underlying algorithms employed inside SuiteSparse:GraphBLAS
// will differ (dot product vs saxpy, for SuiteSparse, for example).
// For example, the Sandia and Sandia2 methods will effectively be
// swapped.

// Method 0 can take a huge amount of memory, for all of A*E.  As a result,
// it often fails for large problems.

// Methods 1 and 2 are much more memory efficient as compare to Method 0,
// taking memory space the same size as A.  But they are slower than methods 3
// to 6.

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

// A is a square symmetric matrix, of any type.  Its values are ignored.  E is
// the edge incidence matrix of A (the values of E are also ignored).
// L=tril(A), and U=triu(A).  See SuiteSparse/GraphBLAS/Demo/tricount.m for a
// complete definition of each method and the matrices A, E, L, and U, and
// citations of relevant references.

// The new GxB_PAIR_INT64 binary operator in SuiteSparse:GraphBLAS v3.2.0 is
// used in the semiring.  This is the function f(x,y)=1, so the values of A, L,
// U, and E are not accessed.  They can have any values and any type.  Only the
// structure of the matrices is used.

// Results are undefined if self-edges exist.

// This code is modified from SuiteSparse/GraphBLAS/Demo/Source/tricount.c.
// It contains no GxB_* extensions and thus it should work in any GraphBLAS
// library.

#define LAGRAPH_FREE_ALL        \
    GrB_free (&S) ;             \
    GrB_free (&C) ;

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// LAGraph_tricount: count the number of triangles in a graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_tricount   // count # of triangles
(
    int64_t *ntri,          // # of triangles
    const int method,       // 0 to 6, see above
    const GrB_Matrix A,     // adjacency matrix for methods 0, 1, and 2
    const GrB_Matrix E,     // edge incidence matrix for method 0
    const GrB_Matrix L,     // L=tril(A) for methods 2, 3, 5, and 6
    const GrB_Matrix U      // U=triu(A) for methods 2, 4, 5, and 6
)
{

    //--------------------------------------------------------------------------
    // check inputs and initialize
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Index n, ne ;
    GrB_Matrix S = NULL, C = NULL ;

#if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
    && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,2,0) )
    // the PAIR function is f(x,y)=1
    GrB_Semiring s = GxB_PLUS_PAIR_INT64 ;
#else
    GrB_Semiring s = LAGraph_PLUS_TIMES_INT64 ;
#endif
    GrB_Monoid sum = LAGraph_PLUS_INT64_MONOID ;
    GrB_Type type = GrB_INT64 ;
    GrB_UnaryOp istwo = LAGraph_ISTWO_INT64 ;

    //--------------------------------------------------------------------------
    // count triangles
    //--------------------------------------------------------------------------

    switch (method)
    {
        case 0:  // minitri:    ntri = nnz (A*E == 2) / 3

            LAGr_Matrix_nrows (&n, A) ;
            LAGr_Matrix_ncols (&ne, E) ;
            LAGr_Matrix_new (&C, type, n, ne) ;
            LAGr_mxm (C, NULL, NULL, s, A, E, NULL) ;
            LAGr_Matrix_new (&S, GrB_BOOL, n, ne) ;
            LAGr_apply (S, NULL, NULL, istwo, C, NULL) ;
            LAGr_reduce (ntri, NULL, sum, S, NULL) ;
            (*ntri) /= 3 ;
            break ;

        case 1:  // Burkhardt:  ntri = sum (sum ((A^2) .* A)) / 6

            LAGr_Matrix_nrows (&n, A) ;
            LAGr_Matrix_new (&C, type, n, n) ;
            LAGr_mxm (C, A, NULL, s, A, A, NULL) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            (*ntri) /= 6 ;
            break ;

        case 2:  // Cohen:      ntri = sum (sum ((L * U) .* A)) / 2

            LAGr_Matrix_nrows (&n, A) ;
            LAGr_Matrix_new (&C, type, n, n) ;
            LAGr_mxm (C, A, NULL, s, L, U, NULL) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            (*ntri) /= 2 ;
            break ;

        case 3:  // Sandia:     ntri = sum (sum ((L * L) .* L))

            // using the masked saxpy3 method
            LAGr_Matrix_nrows (&n, L) ;
            LAGr_Matrix_new (&C, type, n, n) ;
            LAGr_mxm (C, L, NULL, s, L, L, NULL) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        case 4:  // Sandia2:    ntri = sum (sum ((U * U) .* U))

            // using the masked saxpy3 method
            LAGr_Matrix_nrows (&n, U) ;
            LAGr_Matrix_new (&C, type, n, n) ;
            LAGr_mxm (C, U, NULL, s, U, U, NULL) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        case 5:  // SandiaDot:  ntri = sum (sum ((L * U') .* L))

            // using the masked dot product (very efficient)
            LAGr_Matrix_nrows (&n, U) ;
            LAGr_Matrix_new (&C, type, n, n) ;
            LAGr_mxm (C, L, NULL, s, L, U, LAGraph_desc_otoo) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        case 6:  // SandiaDot2: ntri = sum (sum ((U * L') .* U))

            // using the masked dot product (very efficient)
            LAGr_Matrix_nrows (&n, U) ;
            LAGr_Matrix_new (&C, type, n, n) ;
            LAGr_mxm (C, U, NULL, s, U, L, LAGraph_desc_otoo) ;
            LAGr_reduce (ntri, NULL, sum, C, NULL) ;
            break ;

        default:    // invalid method

            return (GrB_INVALID_VALUE) ;
            break ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

