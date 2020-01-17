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
//                             s ≠ i ≠ t
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
// be in their default format (by row); in this case, both phases use a "push"
// direction (a saxpy-based multiply) in SuiteSparse:GraphBLAS.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_WORK                       \
{                                               \
    GrB_free (&frontier) ;                      \
    GrB_free (&paths) ;                         \
    GrB_free (&bc_update) ;                     \
    GrB_free (&W) ;                             \
    if (S != NULL)                              \
    {                                           \
        for (int64_t i = 0 ; i < n ; i++)       \
        {                                       \
            if (S [i] == NULL) break ;          \
            GrB_free (&(S [i])) ;               \
        }                                       \
        free (S) ;                              \
    }                                           \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    GrB_free (centrality) ;         \
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_bc_batch4      // betweeness centrality, batch algorithm
(
    GrB_Vector *centrality,     // centrality(i): betweeness centrality of i
    const GrB_Matrix A,         // input graph, A(i,j) is the edge (i,j)
    const GrB_Matrix AT,        // A'
    const GrB_Index *sources,   // source vertices to compute shortest paths
    int32_t ns                  // number of source vertices
)
{
    GrB_Info info ;
    (*centrality) = NULL ;
    GrB_Index n = 0 ;           // # nodes in the graph

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

    LAGr_Matrix_nrows (&n, A) ;      // # of nodes

    // Create the result vector, one entry for each node
    LAGr_Vector_new (centrality, GrB_FP32, n) ;
    LAGr_Matrix_new (&paths,     GrB_FP32, ns, n) ;
    LAGr_Matrix_new (&frontier,  GrB_FP32, ns, n) ;

    // Initialize paths to source vertices with ones, and the other
    // entries equal to zero.  The paths matrix is dense, and stays that way.
    LAGr_assign (paths, NULL, NULL, 0, GrB_ALL, ns, GrB_ALL, n, NULL) ;
    for (GrB_Index i = 0 ; i < ns ; i++)
    {
        // paths (i,s(i)) = 1
        LAGr_Matrix_setElement (paths, 1, i, sources [i]) ;
        LAGr_Matrix_setElement (frontier, 1, i, sources [i]) ;
    }

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS ) \
        && ( GxB_IMPLEMENTATION >= GxB_VERSION (3,2,0) )
    GrB_Descriptor desc_rc = GrB_DESC_RC ;
    GrB_Descriptor desc_rs = GrB_DESC_RS ;
    GrB_Descriptor desc_t0 = GrB_DESC_T0 ;
    #else
    GrB_Descriptor desc_rc = LAGraph_desc_oocr ;
    GrB_Descriptor desc_rs = LAGraph_desc_ooor ;
    GrB_Descriptor desc_t0 = LAGraph_desc_tooo  ;
    #endif

    // Initial frontier: frontier<!paths>= frontier*A
    LAGr_mxm (frontier, paths, NULL, GxB_PLUS_FIRST_FP32, frontier, A,
        desc_rc) ;

    // Allocate memory for the array of S matrices
    S = (GrB_Matrix *) LAGraph_calloc  (n, sizeof (GrB_Matrix)) ;
    if (S == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    // === Breadth-first search stage ==========================================

    GrB_Index frontier_size ;
    LAGr_Matrix_nvals (&frontier_size, frontier) ;

    int64_t depth ;
    for (depth = 0 ; frontier_size > 0 && depth < n ; depth++)
    {
        // S [depth] = pattern of frontier
        LAGr_Matrix_new (&(S [depth]), GrB_BOOL, ns, n) ;
        LAGr_apply (S [depth], NULL, NULL, GxB_ONE_BOOL, frontier, NULL) ;

        // Accumulate path counts: paths += frontier
        LAGr_assign (paths, NULL, GrB_PLUS_FP32, frontier, GrB_ALL, n, GrB_ALL,
            ns, NULL) ;

        // Update frontier: frontier<!paths> = frontier*A
        LAGr_mxm (frontier, paths, NULL, GxB_PLUS_FIRST_FP32, frontier, A,
            desc_rc) ;

        // Get the size of the current frontier
        LAGr_Matrix_nvals (&frontier_size, frontier) ;
    }

    // === Betweenness centrality computation phase ============================

    // bc_update = ones (ns, n) ; a dense matrix (and stays dense)
    LAGr_Matrix_new (&bc_update, GrB_FP32, ns, n) ;
    LAGr_assign (bc_update, NULL, NULL, 1, GrB_ALL, ns, GrB_ALL, n, NULL) ;
    // W: empty ns-by-n array, as workspace
    LAGr_Matrix_new (&W, GrB_FP32, ns, n) ;

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth-1 ; i > 0 ; i--)
    {
        // Add contributions by successors and mask with that level's frontier

        // W<S[i]> = bc_update ./ path
        LAGr_eWiseMult (W, S [i], NULL, GrB_DIV_FP32, bc_update, paths,
            desc_rs) ;

        // W<S[i−1]> = W * A'
        LAGr_mxm (W, S [i-1], NULL, GxB_PLUS_FIRST_FP32, W, AT, desc_rs) ;

        // bc_update += W .* paths
        // bc_update and paths are both dense, but W is sparse
        LAGr_eWiseMult (bc_update, NULL, GrB_PLUS_FP32, GrB_TIMES_FP32, W,
            paths, NULL) ;
    }

    // Initialize the centrality array with -ns to avoid counting
    // zero length paths
    LAGr_assign (*centrality, NULL, NULL, -ns, GrB_ALL, n, NULL) ;

    // centrality (i) = sum (bc_update (:,i)) for all nodes i
    LAGr_reduce (*centrality, NULL, GrB_PLUS_FP32, GrB_PLUS_FP32, bc_update,
        desc_t0) ;

    LAGRAPH_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

