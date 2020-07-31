//------------------------------------------------------------------------------
// LAGraph_bfs_parent:  push-pull breadth-first search
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

// LAGraph_bfs_parent: direction-optimized push/pull breadth first search,
// contributed by Tim Davis, Texas A&M.  Computes only the BFS tree.
// Requires SuiteSparse:GraphBLAS v4.0.

// Usage:

// info = LAGraph_bfs_parent (&pi, A, AT, source) ;

//      GrB_Vector *pi:  a vector containing the BFS tree, in 1-based indexing.
//          pi(source) = source+1 for source node.  pi(i) = p+1 if p is the
//          parent of i.  If pi is sparse, and pi(i) is not present, then node
//          i has not been reached.  Otherwise, if pi is full, then pi(i)=0
//          indicates that node i was not reached.

//      GrB_Matrix A: a square matrix of any type.  The values of A are not
//          accessed.  The presence of the entry A(i,j) indicates the edge
//          (i,j).  That is, an explicit entry A(i,j)=0 is treated as an edge.

//      GrB_Matrix AT: an optional matrix of any type.  If NULL, the algorithm
//          is a conventional push-only BFS.  If not NULL, AT must be the
//          transpose of A, and a push-pull algorithm is used (NOTE: this
//          assumes GraphBLAS stores its matrix in CSR form; see discussion in
//          LAGraph_bfs_pushpull).  Results are undefined if AT is not NULL but
//          not identical to the transpose of A.

//      int64_t source: the source node for the BFS.

// This algorithm can use the push-pull strategy, which requires both A and
// AT=A' to be passed in.  If the graph is known to be symmetric, then the same
// matrix A can be passed in for both arguments.  Results are undefined if AT
// is not the transpose of A.

// See LAGraph_bfs_pushpull for a discussion of the push/pull strategy.

// References:

// Carl Yang, Aydin Buluc, and John D. Owens. 2018. Implementing Push-Pull
// Efficiently in GraphBLAS. In Proceedings of the 47th International
// Conference on Parallel Processing (ICPP 2018). ACM, New York, NY, USA,
// Article 89, 11 pages. DOI: https://doi.org/10.1145/3225058.3225122

// Scott Beamer, Krste Asanovic and David A. Patterson,
// The GAP Benchmark Suite, http://arxiv.org/abs/1508.03619, 2015.
// http://gap.cs.berkeley.edu/

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&q) ;         \
    GrB_free (&pi) ;        \
}

GrB_Info LAGraph_bfs_parent // push-pull BFS, compute the tree only
(
    // output:
    GrB_Vector *pi_output,  // pi(i) = p+1 if p is the parent of node i
    // inputs:
    GrB_Matrix A,           // input graph, any type
    GrB_Matrix AT,          // transpose of A (optional; push-only if NULL)
    int64_t source          // starting node of the BFS
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION < GxB_VERSION (4,0,0) )
    // SuiteSparse GraphBLAS v4.0 or later required
    return (GrB_INVALID_VALUE) ;
    #else

    GrB_Info info ;
    GrB_Vector q = NULL ;           // nodes visited at each level
    GrB_Vector pi = NULL ;          // parent vector

    if (pi_output == NULL || (A == NULL && AT == NULL))
    {
        // required output argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    bool use_vxm_with_A ;
    GrB_Index nrows, ncols, nvalA ;
    if (A == NULL)
    {
        // only AT is provided
        LAGr_Matrix_ncols (&nrows, AT) ;
        LAGr_Matrix_nrows (&ncols, AT) ;
        LAGr_Matrix_nvals (&nvalA, AT) ;
        use_vxm_with_A = false ;
    }
    else
    {
        // A is provided.  AT may or may not be provided
        LAGr_Matrix_nrows (&nrows, A) ;
        LAGr_Matrix_ncols (&ncols, A) ;
        LAGr_Matrix_nvals (&nvalA, A) ;
        use_vxm_with_A = true ;
    }

    // push/pull requires both A and AT
    bool push_pull = (A != NULL && AT != NULL) ;

    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // check the format of A and AT
    //--------------------------------------------------------------------------

    bool csr = true ;

    // csr is true if A and AT are known (or assumed) to be in CSR format; if
    // false, they are known to be in CSC format.

        GxB_Format_Value A_format = -1, AT_format = -1 ;
        bool A_csr = true, AT_csr = true ;
        if (A != NULL)
        {
            // A_csr is true if accessing A(i,:) is fast
            LAGr_get (A , GxB_FORMAT, &A_format) ;
            A_csr = (A_format == GxB_BY_ROW) ;
        }
        if (AT != NULL)
        {
            // AT_csr is true if accessing AT(i,:) is fast
            LAGr_get (AT, GxB_FORMAT, &AT_format) ;
            AT_csr = (AT_format == GxB_BY_ROW) ;
        }
        // Assume CSR if A(i,:) and AT(i,:) are both fast.  If csr is false,
        // then the algorithm below will reverse the use of vxm and mxv.
        csr = A_csr && AT_csr ;
        if (push_pull)
        {
            // both A and AT are provided.  Require they have the same format.
            // Either both A(i,:) and AT(i,:) are efficient to accesss, or both
            // A(:,j) and AT(:,j) are efficient to access.
            if (A_csr != AT_csr)
            {
                LAGRAPH_ERROR ("A and AT must in the same format:\n"
                    "both GxB_BY_ROW, or both GxB_BY_COL",
                    GrB_INVALID_VALUE) ;
            }
        }
#if 0
        else
        {
            // only A or AT are provided.  Refuse to do the pull-only version.
            if (A != NULL && A_format == GxB_BY_COL)
            {
                // this would result in a pull-only BFS ... exceedingly slow
                LAGRAPH_ERROR (
                    "SuiteSparse: AT not provided, so A must be GxB_BY_ROW\n"
                    "(or provide both A and AT, both in the same format,\n"
                    "either both GxB_BY_COL or both GxB_BY_ROW)",
                    GrB_INVALID_VALUE) ;
            }
            if (AT != NULL && AT_format == GxB_BY_ROW)
            {
                // this would result in a pull-only BFS ... exceedingly slow
                LAGRAPH_ERROR (
                    "SuiteSparse: A not provided, so AT must be GxB_BY_COL\n"
                    "(or provide both A and AT, both in the same format,\n"
                    "either both GxB_BY_COL or both GxB_BY_ROW)",
                    GrB_INVALID_VALUE) ;
            }
        }
#endif

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n = nrows ;
    GrB_Type int_type = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    GrB_Semiring semiring ;

    // create an sparse integer vector q, and set q(source) to source+1
    LAGr_Vector_new (&q, int_type, n) ;
    LAGr_Vector_setElement (q, source+1, source) ;
    GrB_Index nq = 1 ;          // number of nodes in the current level

    if (n > INT32_MAX)
    {
        semiring = GxB_ANY_SECONDI1_INT64 ;
    }
    else
    {
        semiring = GxB_ANY_SECONDI1_INT32 ;
    }

    // pi = a dense vector of all zeros
    LAGr_Vector_new (&pi, int_type, n) ;
    LAGr_assign (pi, NULL, NULL, 0, GrB_ALL, n, NULL) ;

    // pi (source) = source+1 denotes a root of the BFS tree
    LAGr_Vector_setElement (pi, source+1, source) ;

    // average node degree
    double d = (n == 0) ? 0 : (((double) nvalA) / (double) n) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    for (int64_t nvisited = 0 ; nvisited < n ; nvisited += nq)
    {

        //----------------------------------------------------------------------
        // select push vs pull
        //----------------------------------------------------------------------

        if (push_pull)
        {
            double pushwork = d * nq ;
            double expected = (double) n / (double) (nvisited+1) ;
            double per_dot = LAGRAPH_MIN (d, expected) ;
            double binarysearch = (3 * (1 + log2 ((double) nq))) ;
            double pullwork = (n-nvisited) * per_dot * binarysearch ;
            use_vxm_with_A = (pushwork < pullwork) ;
            if (!csr)
            {
                // Neither A(i,:) nor AT(i,:) is efficient.  Instead, both
                // A(:,j) and AT(:,j) is fast (that is, the two matrices
                // are in CSC format).  Swap vxm and mxv.
                use_vxm_with_A = !use_vxm_with_A ;
            }
        }

        //----------------------------------------------------------------------
        // q = next level of the BFS
        //----------------------------------------------------------------------

double t = omp_get_wtime ( ) ;

        if (use_vxm_with_A)
        {
            // q'<!pi> = q'*A
            // this is a push step if A is in CSR format; pull if CSC
            LAGr_vxm (q, pi, NULL, semiring, q, A, GrB_DESC_RC) ;
        }
        else
        {
            // q<!pi> = AT*q
            // this is a pull step if AT is in CSR format; push if CSC
            LAGr_mxv (q, pi, NULL, semiring, AT, q, GrB_DESC_RC) ;
        }

t = omp_get_wtime ( ) - t ;

        LAGr_Vector_nvals (&nq, q) ;
printf ("%s %10ld %12.5f\n", use_vxm_with_A ? "push" : "pull", nq, t) ;
        if (nq == 0) break ;

        //----------------------------------------------------------------------
        // assign parents
        //----------------------------------------------------------------------

        // q(i) currently contains the parent+1 of node i in tree (off by one
        // so it won't have any zero values, for valued mask).
        // pi<q> = q
        LAGr_assign (pi, q, NULL, q, GrB_ALL, n, GrB_DESC_S) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*pi_output) = pi ;
    pi = NULL ;
    LAGRAPH_FREE_ALL ;      // free all workspace (except for result pi)
    return (GrB_SUCCESS) ;
    #endif
}

