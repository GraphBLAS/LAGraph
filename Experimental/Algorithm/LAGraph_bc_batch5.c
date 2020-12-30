//------------------------------------------------------------------------------
// LAGraph_bc_batch: Brandes' algorithm for computing betweeness centrality
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

// LAGraph_bc_batch: Batch algorithm for computing betweeness centrality.
// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University.
// Adapted from GraphBLAS C API Spec, Appendix B.4.

// LAGraph_bc_batch computes an approximation of the betweenness centrality of
// all nodes in a graph using a batched version of Brandes' algorithm.
//                               ____
//                               \      sigma(s,t | i)
//    Betweenness centrality =    \    ----------------
//           of node i            /       sigma(s,t)
//                               /___
//                            s != i != t
//
// Where sigma(s,t) is the total number of shortest paths from node s to
// node t, and sigma(s,t | i) is the total number of shortest paths from
// node s to node t that pass through node i.
//
// Note that the true betweenness centrality requires computing shortest paths
// from all nodes s to all nodes t (or all-pairs shortest paths), which can be
// expensive to compute. By using a reasonably sized subset of source nodes, an
// approximation can be made.
//
// LAGraph_bc_batch performs simultaneous breadth-first searches of the entire
// graph starting at a given set of source nodes. This pass discovers all
// shortest paths from the source nodes to all other nodes in the graph. After
// the BFS is complete, the number of shortest paths that pass through a given
// node is tallied by reversing the traversal. From this, the (approximate)
// betweenness centrality is computed.

// A represents the graph, and AT must equal A'.  A must be square, and can be
// unsymmetric.  Self-edges are OK.  The values of A and AT are ignored; just
// the pattern of two matrices are used.  For best performance, A and AT should
// be in their default format (by row).

// Each phases uses push-pull direction optimization.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_WORK                       \
{                                               \
    LAGr_free (&frontier) ;                     \
    LAGr_free (&paths) ;                        \
    LAGr_free (&bc_update) ;                    \
    LAGr_free (&W) ;                            \
    if (S != NULL)                              \
    {                                           \
        for (int64_t i = 0 ; i < n ; i++)       \
        {                                       \
            if (S [i] == NULL) break ;          \
            LAGr_free (&(S [i])) ;              \
        }                                       \
        free (S) ;                              \
    }                                           \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    LAGr_free (centrality) ;        \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_bc_batch5      // betweeness centrality, batch algorithm
(
    GrB_Vector *centrality,     // centrality(i): betweeness centrality of i
    const GrB_Matrix A,         // input graph, A(i,j) is the edge (i,j)
    const GrB_Matrix AT,        // A'
    const GrB_Index *sources,   // source vertices to compute shortest paths
    int32_t ns                  // number of source vertices
)
{

    // =========================================================================
    // === initializations =====================================================
    // =========================================================================

    double tt1 = omp_get_wtime ( ) ;
    GrB_Info info ;
    (*centrality) = NULL ;
    GrB_Index n ;                   // # nodes in the graph

    // Array of BFS search matrices.
    // S [i] is a sparse matrix that stores the depth at which each vertex is
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Matrix *S = NULL ;

    // Frontier matrix, a sparse matrix.
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Matrix frontier = NULL ;

    // Paths matrix holds the number of shortest paths for each node and
    // starting node discovered so far.  A dense matrix that is updated with
    // sparse updates, and also used as a mask.
    GrB_Matrix paths = NULL ;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node.  A dense matrix.
    GrB_Matrix bc_update = NULL ;

    // Temporary workspace matrix (sparse).
    GrB_Matrix W = NULL ;

    // Initialize paths and frontier with source notes
    LAGr_Matrix_nrows (&n, A) ;
    LAGr_Matrix_new (&paths,     GrB_FP32, ns, n) ;
    LAGr_Matrix_new (&frontier,  GrB_FP32, ns, n) ;
    GxB_set (paths, GxB_SPARSITY_CONTROL, GxB_BITMAP + GxB_FULL) ;
    for (GrB_Index i = 0 ; i < ns ; i++)
    {
        // paths (i,s(i)) = 1
        LAGr_Matrix_setElement (paths,    1, i, sources [i]) ;
        LAGr_Matrix_setElement (frontier, 1, i, sources [i]) ;
    }

    // Initial frontier: frontier<!paths>= frontier*A
    LAGr_mxm (frontier, paths, NULL, GxB_PLUS_FIRST_FP32, frontier, A,
        GrB_DESC_RC) ;

    // Allocate memory for the array of S matrices
    S = (GrB_Matrix *) LAGraph_calloc  (n, sizeof (GrB_Matrix)) ;
    if (S == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    // =========================================================================
    // === Breadth-first search stage ==========================================
    // =========================================================================

    bool last_was_pull = false ;
    GrB_Index frontier_size, last_frontier_size = 0 ;
    LAGr_Matrix_nvals (&frontier_size, frontier) ;

    tt1 = omp_get_wtime ( ) - tt1 ;
    printf ("init phase:    %.3g sec\n", tt1) ;
    tt1 = omp_get_wtime ( ) ;

    int64_t depth ;
    for (depth = 0 ; frontier_size > 0 && depth < n ; depth++)
    {

        //----------------------------------------------------------------------
        // S [depth] = pattern of frontier
        //----------------------------------------------------------------------

//      double ttt = omp_get_wtime ( ) ;
        LAGr_Matrix_new (&(S [depth]), GrB_BOOL, ns, n) ;
        LAGr_apply (S [depth], NULL, NULL, GxB_ONE_BOOL, frontier, NULL) ;

        // TODO (in SuiteSparse:GraphBLAS, not LAGraph): constructing S will be
        // faster with uniform-valued matrices, once they are added to
        // SuiteSparse:GraphBLAS.

        // TODO (in SuiteSparse:GraphBLAS, not LAGraph): mxm can produce a
        // jumbled matrix (frontier), which is then sorted by the apply (here)
        // and the eMult (below).  Both the apply and the eMult could tolerate
        // jumbled inputs, and produce jumbled outputs, so if enough
        // jumbled-exploit is added, the lazy sort could be so lazy that it
        // never happens (like the BFS).

        //----------------------------------------------------------------------
        // Accumulate path counts: paths += frontier
        //----------------------------------------------------------------------

        LAGr_assign (paths, NULL, GrB_PLUS_FP32, frontier, GrB_ALL, n, GrB_ALL,
            ns, NULL) ;

        //----------------------------------------------------------------------
        // Update frontier: frontier<!paths> = frontier*A
        //----------------------------------------------------------------------

        double frontier_density = ((double) frontier_size) / (double) (ns*n) ;
        bool growing = frontier_size > last_frontier_size ;
        // pull if frontier is more than 10% dense,
        // or > 6% dense and last step was pull
        bool do_pull = (frontier_density > 0.10) ||
                      ((frontier_density > 0.06) && last_was_pull) ;
//      // pull if frontier is more than 4% dense and growing, or
//      // more than 10% dense yet shrinking.
//      bool do_pull = frontier_density > ((growing) ? 0.04 : 0.10) ;
//      // pull if frontier is more than 10% dense
//      bool do_pull = (frontier_density > 0.10) ;

        if (do_pull)
        {
            GxB_set (frontier, GxB_SPARSITY_CONTROL, GxB_BITMAP) ;
            LAGr_mxm (frontier, paths, NULL, GxB_PLUS_FIRST_FP32, frontier, AT,
                GrB_DESC_RCT1) ;
        }
        else // push
        {
            GxB_set (frontier, GxB_SPARSITY_CONTROL, GxB_SPARSE) ;
            LAGr_mxm (frontier, paths, NULL, GxB_PLUS_FIRST_FP32, frontier, A,
                GrB_DESC_RC) ;
        }

        //----------------------------------------------------------------------
        // Get the size of the current frontier
        //----------------------------------------------------------------------

        last_frontier_size = frontier_size ;
        last_was_pull = do_pull ;
        LAGr_Matrix_nvals (&frontier_size, frontier) ;

//      ttt = omp_get_wtime ( ) - ttt ;
//      printf ("frontier_size %12ld %8.2f%% growing: %d   %s %0.3g sec\n",
//          last_frontier_size, 100 * frontier_density, growing,
//          do_pull ? "pull" : "    ", ttt) ;
    }

    LAGr_free (&frontier) ;

    tt1 = omp_get_wtime ( ) - tt1 ;
    printf ("forward phase: %.3g sec\n", tt1) ;
    tt1 = omp_get_wtime ( ) ;

    // =========================================================================
    // === Betweenness centrality computation phase ============================
    // =========================================================================

    // bc_update = ones (ns, n) ; a full matrix (and stays full)
    LAGr_Matrix_new (&bc_update, GrB_FP32, ns, n) ;
    LAGr_assign (bc_update, NULL, NULL, 1, GrB_ALL, ns, GrB_ALL, n, NULL) ;
    // W: empty ns-by-n array, as workspace
    LAGr_Matrix_new (&W, GrB_FP32, ns, n) ;

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth-1 ; i > 0 ; i--)
    {
//      double ttt = omp_get_wtime ( ) ;

        //----------------------------------------------------------------------
        // W<S[i]> = bc_update ./ paths
        //----------------------------------------------------------------------

        // Add contributions by successors and mask with that level's frontier
        LAGr_eWiseMult (W, S [i], NULL, GrB_DIV_FP32, bc_update, paths,
            GrB_DESC_RS) ;

        //----------------------------------------------------------------------
        // W<S[i−1]> = W * A'
        //----------------------------------------------------------------------

        // pull if W is more than 10% dense and nnz(W)/nnz(S[i-1]) > 1
        // or if W is more than 1% dense and nnz(W)/nnz(S[i-1]) > 10
        GrB_Index wsize, ssize ;
        GrB_Matrix_nvals (&wsize, W) ;
        GrB_Matrix_nvals (&ssize, S [i-1]) ;
        double w_density    = ((double) wsize) / ((double) (ns*n)) ;
        double w_to_s_ratio = ((double) wsize) / ((double) ssize) ;
        bool do_pull = (w_density > 0.1  && w_to_s_ratio > 1.) ||
                       (w_density > 0.01 && w_to_s_ratio > 10.) ;

        if (do_pull)
        {
            GxB_set (W, GxB_SPARSITY_CONTROL, GxB_BITMAP) ;
            LAGr_mxm (W, S [i-1], NULL, GxB_PLUS_FIRST_FP32, W, A,
                GrB_DESC_RST1) ;
        }
        else // push
        {
            GxB_set (W, GxB_SPARSITY_CONTROL, GxB_SPARSE) ;
            LAGr_mxm (W, S [i-1], NULL, GxB_PLUS_FIRST_FP32, W, AT,
                GrB_DESC_RS) ;
        }

        //----------------------------------------------------------------------
        // bc_update += W .* paths
        //----------------------------------------------------------------------

        // bc_update is full, paths is bitmap/full, W is sparse/bitmap
        LAGr_eWiseMult (bc_update, NULL, GrB_PLUS_FP32, GrB_TIMES_FP32, W,
            paths, NULL) ;

//      ttt = omp_get_wtime ( ) - ttt ;
//      printf ("back wsize    %12ld %8.2f%% w/s %8.2f %s %0.3g sec\n",
//          wsize, 100 * w_density, w_to_s_ratio,
//          do_pull ? "pull" : "    ", ttt) ;
    }

    // =========================================================================
    // === finalize the centrality =============================================
    // =========================================================================

    // Initialize the centrality array with -ns to avoid counting
    // zero length paths
    LAGr_Vector_new (centrality, GrB_FP32, n) ;
    LAGr_assign (*centrality, NULL, NULL, -ns, GrB_ALL, n, NULL) ;

    // centrality (i) = sum (bc_update (:,i)) for all nodes i
    LAGr_reduce (*centrality, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32, bc_update,
        GrB_DESC_T0) ;

    LAGRAPH_FREE_WORK ;
    tt1 = omp_get_wtime ( ) - tt1 ;
    printf ("back phase:    %.3g sec\n", tt1) ;
    return (GrB_SUCCESS) ;
}

