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

#include "../../Source/Utility/LAGraph_internal.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&w) ;         \
    GrB_free (&q) ;         \
    GrB_free (&pi) ;        \
}

GrB_Info bfs_log_parent // push-pull BFS, compute the tree only
(
    // output:
    GrB_Vector *pi_output,  // pi(i) = p+1 if p is the parent of node i
    // inputs:
    GrB_Matrix A,           // input graph, any type
    GrB_Matrix AT,          // transpose of A (optional; push-only if NULL)
    GrB_Vector Degree,      // Degree(i) is the out-degree of node i
                            // (optional: push-only if NULL)
    int64_t source          // starting node of the BFS
    , FILE *file
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
    GrB_Vector q = NULL ;           // the current frontier
    GrB_Vector pi = NULL ;          // parent vector
    GrB_Vector w = NULL ;           // to compute work remaining

    if (pi_output == NULL || (A == NULL && AT == NULL))
    {
        // required output argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    GrB_Index nrows, ncols, nvals ;
    if (A == NULL)
    {
        // only AT is provided
        LAGr_Matrix_ncols (&nrows, AT) ;
        LAGr_Matrix_nrows (&ncols, AT) ;
        LAGr_Matrix_nvals (&nvals, AT) ;
    }
    else
    {
        // A is provided.  AT may or may not be provided
        LAGr_Matrix_nrows (&nrows, A) ;
        LAGr_Matrix_ncols (&ncols, A) ;
        LAGr_Matrix_nvals (&nvals, A) ;
    }

    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE) ;
    }

    //--------------------------------------------------------------------------
    // check the format of A and AT
    //--------------------------------------------------------------------------

    bool A_csr = true, AT_csr = true ;
    if (A != NULL)
    {
        GxB_Format_Value A_format ;
        LAGr_get (A , GxB_FORMAT, &A_format) ;
        A_csr = (A_format == GxB_BY_ROW) ;
    }
    if (AT != NULL)
    {
        GxB_Format_Value AT_format ;
        LAGr_get (AT, GxB_FORMAT, &AT_format) ;
        AT_csr = (AT_format == GxB_BY_ROW) ;
    }

    bool vxm_is_push = (A  != NULL &&  A_csr ) ;    // vxm (q,A) is a push step
    bool vxm_is_pull = (A  != NULL && !A_csr ) ;    // vxm (q,A) is a pull step

    bool mxv_is_push = (AT != NULL && !AT_csr) ;    // mxv (AT,q) is a push step
    bool mxv_is_pull = (AT != NULL &&  AT_csr) ;    // mxv (AT,q) is a pull step

    // can_push is true if the push-step can be performed
    bool can_push = vxm_is_push || mxv_is_push ;

    // can_pull is true if the pull-step can be performed
    bool can_pull = vxm_is_pull || mxv_is_pull ;

    // direction-optimization requires both push and pull, and Degree
    bool push_pull = can_push && can_pull && (Degree != NULL) ;

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

    if (push_pull) LAGr_Vector_new (&w, GrB_INT64, n) ;

fprintf (file, "\nk = k+1 ; s{k} = %lu ;\n", source) ;
fprintf (file, "results{k} = [\n") ;

    double alpha = 4.0 ;
    double beta1 = 8.0 ;
    double beta2 = 500.0 ;
    int64_t n_over_beta1 = (int64_t) (((double) n) / beta1) ;
    int64_t n_over_beta2 = (int64_t) (((double) n) / beta2) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    bool do_push = can_push ;   // start with push, if available
    GrB_Index last_nq = 0 ;
    int64_t edges_unexplored = nvals ;
    bool any_pull = false ;     // true if any pull phase has been done

    for (int64_t nvisited = 0 ; nvisited < n ; nvisited += nq)
    {

        //----------------------------------------------------------------------
        // select push vs pull
        //----------------------------------------------------------------------

        int64_t edges_in_frontier = 0 ;
        if (push_pull)
        {
            if (do_push && can_pull)
            {
                // check for switch from push to pull
                bool growing = nq > last_nq ;
                bool switch_to_pull ;
                if (any_pull)
                { 
                    // once any pull phase has been done, the # of edges in the
                    // frontier has no longer been tracked.  But now the BFS
                    // has switched back to push, and we're checking for yet
                    // another switch to pull.  This switch is unlikely, so
                    // just keep track of the size of the frontier, and switch
                    // if it starts growing again and is getting big.
                    switch_to_pull = (growing && nq > n_over_beta1) ;
                }
                else
                {
                    // update the # of unexplored edges
                    // w<q>=Degree
                    // w(i) = outdegree of node i if node i is in the queue
                    // GxB_set ((GrB_Matrix) q, GxB_SPARSITY, GxB_SPARSE) ;
                    LAGr_assign (w, q, NULL, Degree, GrB_ALL, n, GrB_DESC_RS) ;
                    // edges_in_frontier = sum (w) = # of edges incident on all
                    // nodes in the current frontier
                    LAGr_reduce (&edges_in_frontier, NULL,
                        GrB_PLUS_MONOID_INT64, w, NULL) ;
                    edges_unexplored -= edges_in_frontier ;
                    switch_to_pull = growing &&
                        (edges_in_frontier > (edges_unexplored / alpha)) ;
                }
                if (switch_to_pull)
                {
                    // the # of edges incident on 
                    do_push = false ;
                }
            }
            else if (!do_push && can_push)
            {
                // check for switch from pull to push
                bool shrinking = nq < last_nq ;
                if (shrinking && (nq < n_over_beta2))
                {
                    do_push = true ;
                }
            }
        }
        any_pull = any_pull || (!do_push) ;

        //----------------------------------------------------------------------
        // q = next level of the BFS
        //----------------------------------------------------------------------

        GxB_set ((GrB_Matrix) q, GxB_SPARSITY,  
            do_push ? GxB_SPARSE : GxB_BITMAP) ;

double t ;

        if ((do_push && vxm_is_push) || (!do_push && vxm_is_pull))
        {
            // q'<!pi> = q'*A
            // this is a push if A is in CSR format; pull if A is in CSC
t = omp_get_wtime ( ) ;
            LAGr_vxm (q, pi, NULL, semiring, q, A, GrB_DESC_RC) ;
        }
        else // ((!do_push && mxv_is_pull) || (do_push && mxv_is_push))
        {
            // q<!pi> = AT*q
            // this is a pull if AT is in CSR format; push if AT is in CSC
t = omp_get_wtime ( ) ;
            LAGr_mxv (q, pi, NULL, semiring, AT, q, GrB_DESC_RC) ;
        }

t = omp_get_wtime ( ) - t ;
fprintf (file, "%d  %16lu %16lu %16ld    %g\n",
    do_push, nq, nvisited, edges_in_frontier, t) ;


        last_nq = nq ;
        LAGr_Vector_nvals (&nq, q) ;
        if (nq == 0) break ;

        //----------------------------------------------------------------------
        // assign parents
        //----------------------------------------------------------------------

        // q(i) currently contains the parent+1 of node i in tree.
        // pi<q> = q
        LAGr_assign (pi, q, NULL, q, GrB_ALL, n, GrB_DESC_S) ;
    }

fprintf (file, "] ;\n") ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*pi_output) = pi ;
    pi = NULL ;
    LAGRAPH_FREE_ALL ;      // free all workspace (except for result pi)
    return (GrB_SUCCESS) ;
    #endif
}

