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

// A_matrix represents the graph.  It must be square, and can be unsymmetric.
// Self-edges are OK.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_WORK                       \
{                                               \
    GrB_free(&frontier);                        \
    GrB_free(&paths);                           \
    GrB_free(&bc_update);                       \
    GrB_free(&temp);                            \
    if (S_array != NULL)                        \
    {                                           \
        for (int64_t i = 0; i < n ; i++)        \
        {                                       \
            if (S_array [i] == NULL) break ;    \
            GrB_free (&(S_array [i])) ;         \
        }                                       \
        free (S_array) ;                        \
    }                                           \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    GrB_free (centrality) ;         \
}

#if 0
// select FP64
#define REAL_t                  double
#define LAGr_REAL_TYPE          GrB_FP64
#define LAGr_PLUS_SECOND_REAL   GxB_PLUS_SECOND_FP64
#define LAGr_PLUS_REAL          GrB_PLUS_FP64
#define LAGr_TIMES_REAL         GrB_TIMES_FP64
#define LAGr_DIV_REAL           GrB_DIV_FP64
#else
// select FP32
#define REAL_t                  float
#define LAGr_REAL_TYPE          GrB_FP32
#define LAGr_PLUS_SECOND_REAL   GxB_PLUS_SECOND_FP32
#define LAGr_PLUS_REAL          GrB_PLUS_FP32
#define LAGr_TIMES_REAL         GrB_TIMES_FP32
#define LAGr_DIV_REAL           GrB_DIV_FP32
#endif


GrB_Info LAGraph_bc_batch3 // betweeness centrality, batch algorithm
(
    GrB_Vector *centrality,    // centrality(i) is the betweeness centrality of node i
    const GrB_Matrix A_matrix, // input graph, treated as if boolean in semiring
    const GrB_Matrix AT_matrix, // A'
    const GrB_Index *sources,  // source vertices from which to compute shortest paths
    int32_t num_sources        // number of source vertices (length of s)
)
{
    GrB_Info info ;
    (*centrality) = NULL ;
    GrB_Index n = 0 ; // Number of nodes in the graph

    // Array of BFS search matrices
    // S_array[i] is a matrix that stores the depth at which each vertex is 
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Matrix *S_array = NULL ;

    // Frontier matrix
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Matrix frontier = NULL ;
    
    // Paths matrix holds the number of shortest paths for each node and 
    // starting node discovered so far.
    GrB_Matrix paths = NULL ;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node
    GrB_Matrix bc_update = NULL ;

    // Temporary workspace matrix
    GrB_Matrix temp = NULL ;

    GxB_Format_Value a_fmt, at_fmt ;
    GxB_get (A_matrix,  GxB_FORMAT, &a_fmt ) ;
    GxB_get (AT_matrix, GxB_FORMAT, &at_fmt) ;
    if (a_fmt != GxB_BY_ROW || at_fmt != GxB_BY_ROW)
    {
        LAGRAPH_ERROR ("A and AT must be stored by row", GrB_INVALID_VALUE) ;
    }

    int64_t depth = 0; // Initial BFS depth
    LAGr_Matrix_nrows(&n, A_matrix); // Get dimensions

    // Create the result vector, one entry for each node
    LAGr_Vector_new(centrality, LAGr_REAL_TYPE, n);

    // Initialize paths to source vertices with ones
    // paths[s[i],i]=1 for i=[0, ..., num_sources)
    LAGr_Matrix_new(&paths, LAGr_REAL_TYPE, n, num_sources);
    GxB_set (paths, GxB_FORMAT, GxB_BY_COL) ;
    for (GrB_Index i = 0; i < num_sources; ++i)
    {
        // paths [s[i],i] = 1
        LAGr_Matrix_setElement(paths, (int64_t) 1, sources[i], i);
    }

    // Create frontier matrix and initialize to outgoing nodes from
    // all source nodes
    LAGr_Matrix_new(&frontier, LAGr_REAL_TYPE, n, num_sources);
    GxB_set (frontier, GxB_FORMAT, GxB_BY_COL) ;
    for (GrB_Index i = 0; i < num_sources; ++i)
    {
        // frontier [s[i],i] = 1 ;
        LAGr_Matrix_setElement(frontier, (int64_t) 1, sources[i], i);
    }

    // AT = A'
    // frontier <!paths> = AT (:,sources)
    // LAGr_extract(frontier, paths, NULL, A_matrix, GrB_ALL, n, sources, num_sources, LAGraph_desc_tocr);

    // Initial frontier: frontier<!paths>=A’ +.∗ frontier
    LAGr_mxm(frontier, paths, NULL, LAGr_PLUS_SECOND_REAL, A_matrix, frontier, LAGraph_desc_tocr);

    // Allocate memory for the array of S matrices
    S_array = (GrB_Matrix*) calloc(n, sizeof(GrB_Matrix));
    if (S_array == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    GxB_print (paths, 3) ;
    GxB_print (frontier, 3) ;

    // === Breadth-first search stage ==========================================
    GrB_Index frontier_size = 0;
    do
    {
        printf ("depth::: %g\n", (double) depth) ;

        // Create the current search matrix - one column for each source/BFS
        LAGr_Matrix_new(&(S_array[depth]), GrB_BOOL, n, num_sources);
        GxB_set (S_array [depth], GxB_FORMAT, GxB_BY_COL) ;

        // Copy the current frontier to S
        LAGr_apply(S_array[depth], NULL, NULL, GrB_IDENTITY_BOOL, frontier, NULL);

        // Accumulate path counts: paths += frontier
        LAGr_assign(paths, NULL, LAGr_PLUS_REAL, frontier, GrB_ALL, n, GrB_ALL, num_sources, NULL);

        // Update frontier: frontier<!paths>=A’ +.∗ frontier
        LAGr_mxm(frontier, paths, NULL, LAGr_PLUS_SECOND_REAL, A_matrix, frontier, LAGraph_desc_tocr);

        // Get the size of the current frontier
        LAGr_Matrix_nvals(&frontier_size, frontier);

        GxB_print (paths, 3) ;
        GxB_print (frontier, 3) ;

        depth = depth + 1;

    } while (frontier_size > 0); // Repeat frontier is empty

    // === Betweenness centrality computation phase ============================

    // Create the update matrix and initialize it to 1
    LAGr_Matrix_new(&bc_update, LAGr_REAL_TYPE, n, num_sources);
    GxB_set (bc_update, GxB_FORMAT, GxB_BY_COL) ;
    LAGr_assign(bc_update, NULL, NULL, 1.0f, GrB_ALL, n, GrB_ALL, num_sources, NULL);

    LAGr_Matrix_new(&temp, LAGr_REAL_TYPE, n, num_sources);
    GxB_set (temp, GxB_FORMAT, GxB_BY_COL) ;

    printf ("\n ---------------------- backtrack\n") ;

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth - 1; i > 0; i--)
    {
        // Add contributions by successors and mask with that BFS level’s frontier
        printf ("\nback depth %g\n", (double) i) ;

        // temp<S_array[i]> = (1 ./ nsp) .∗ bc_update
        LAGr_eWiseMult(temp, S_array[i], NULL, LAGr_DIV_REAL, bc_update, paths, LAGraph_desc_ooor);
        GxB_print (temp, 3) ;

        // temp<S_array[i−1]> = (AT' ∗ temp), to use saxpy
        LAGr_mxm (temp, S_array[i-1], NULL, LAGr_PLUS_SECOND_REAL, AT_matrix, temp, LAGraph_desc_toor);
        GxB_print (temp, 3) ;

        // bc_update += temp .∗ paths
        LAGr_eWiseMult(bc_update, NULL, LAGr_PLUS_REAL, LAGr_TIMES_REAL, temp, paths, NULL);
        GxB_print (bc_update, 3) ;

    }
    // Initialize the centrality array with -(num_sources) to avoid counting
    // zero length paths
    LAGr_assign(*centrality, NULL, NULL, -(REAL_t)num_sources, GrB_ALL, n, NULL);

    // centrality += update
    LAGr_reduce(*centrality, NULL, LAGr_PLUS_REAL, LAGr_PLUS_REAL, bc_update, NULL);

    LAGRAPH_FREE_WORK ;
    return GrB_SUCCESS;
}
