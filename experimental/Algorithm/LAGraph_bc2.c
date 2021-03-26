//------------------------------------------------------------------------------
// LAGraph_bc: Brandes' algorithm for computing betweeness centrality
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

// LAGraph_bc: Brandes' algorithm for computing betweeness centrality.
// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University.
// Adapted from GraphBLAS C API Spec, Appendix B.3.

// LAGraph_bc computes an approximation of the betweenness centrality of all 
// nodes in a graph using Brandes' algorithm. 
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
// expensive to compute. By using a single source node, an approximation can be
// made, and by repeatedly calling LAGraph_bc from a set of different source
// nodes, a progressively better approximation can be obtained.
//
// LAGraph_bc performs a breadth-first search of the entire graph starting at
// a given source node. This pass discovers all shortest paths from the source
// node to all other nodes in the graph. After the BFS is complete, the number
// of shortest paths that pass through a given node is tallied by reversing the
// traversal. From this, the (approximate) betweenness centrality is computed.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL                      \
{                                             \
    GrB_free (&S_matrix);                     \
    GrB_free (&frontier);                     \
    GrB_free (&paths);                        \
    GrB_free (&structural_complement_output); \
    GrB_free (&transpose_first_arg);          \
    GrB_free (&temp1);                        \
    GrB_free (&temp2);                        \
}

GrB_Info LAGraph_bc2     // betweeness centrality
(
    GrB_Vector *centrality, // centrality(i): betweeness centrality of node i
    GrB_Matrix A_matrix,    // input graph
    GrB_Index source        // source vertex
)
{

    GrB_Index n; // Number of nodes in the graph
    
    // BFS search matrix
    // Stores depth at which each vertex is seen
    GrB_Matrix S_matrix;

    // Frontier vector
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Vector frontier;

    // Vector storing number of shortest paths to each vertex
    // Updated every BFS iteration as more paths are explored
    GrB_Vector paths;

    GrB_Descriptor structural_complement_output;
    GrB_Descriptor transpose_first_arg;
    
    // Temp vectors used in computing the centrality update
    GrB_Vector temp1, temp2;

    LAGr_Matrix_nrows(&n, A_matrix); // Get dimensions

    // Create the result vector for storing the final centrality metric
    LAGr_Vector_new(centrality, GrB_FP64, n);
    
    // Create the search matrix. Dimensions are n x n:
    //  n vertices and (worst case) n levels in the BFS
    LAGr_Matrix_new(&S_matrix, GrB_FP64, n, n);
    
    // Create the frontier vector - one entry for each vertex
    LAGr_Vector_new(&frontier, GrB_FP64, n);

    // Initialize the frontier vector with the source vertex
    LAGr_Vector_setElement(frontier, 1, source);
    
    // Initialize the number of paths with a copy of the frontier vector
    LAGr_Vector_dup(&paths, frontier);
    
    // Create a descriptor to return the structural complement of the output
    LAGr_Descriptor_new(&structural_complement_output);
    LAGr_Descriptor_set(structural_complement_output, GrB_MASK, GrB_SCMP);
    LAGr_Descriptor_set(structural_complement_output, GrB_OUTP, GrB_REPLACE);
    
    // Create a descriptor to use the transpose of the first input argument
    LAGr_Descriptor_new(&transpose_first_arg);
    LAGr_Descriptor_set(transpose_first_arg, GrB_INP0, GrB_TRAN);

    // === Breadth-first search stage ==========================================
    int64_t depth = 0;  // Start at depth 0
    double sum = 0;    // Sum of shortest paths to vertices at current depth
                        // Equal to sum(frontier). Continue BFS until new paths
                        //  are no shorter than any existing paths.
    do
    {
        // Set this column of the S_matrix matrix equal to the frontier vector
        // S_matrix(:,depth) = frontier
        LAGr_assign(S_matrix, GrB_NULL, GrB_NULL, frontier, depth, GrB_ALL, n, GrB_NULL);

        // Traverse to the next level of the BFS.
        // frontier = frontier*S_matrix masked by paths
        // The structural complement descriptor is used to denote !paths
        LAGr_vxm(frontier, paths, GrB_NULL, GxB_PLUS_TIMES_FP64, frontier, A_matrix, structural_complement_output);

        // Accumulate shortest paths: paths = paths + frontier
        LAGr_eWiseAdd(paths, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_FP64, paths, frontier, GrB_NULL);

        // Sum path counts: sum(frontier)
        LAGr_reduce(&sum, GrB_NULL, GxB_PLUS_FP64_MONOID, frontier, GrB_NULL);

        // Increase BFS depth by one
        depth = depth + 1;

    } while (sum > 0); // Repeat until sum(frontier) = 0

    // === Betweenness centrality computation phase ============================
    LAGr_Vector_new(&temp1, GrB_FP64, n);
    LAGr_Vector_new(&temp2, GrB_FP64, n);

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth; i > 1; i--)
    {
        // TODO: Optimize this - probably possible with only one temp vector
        
        // Update formula for the centrality update is shown here:
        // u = S_matrix(i-1,:)' .* A_matrix * (1 + centrality) ./ S_matrix(i,:)'
        // We need to build this complex operation piecewise.

        // temp1 = ones(1,n)
        LAGr_assign(temp1, GrB_NULL, GrB_NULL, 1.0f, GrB_ALL, n, GrB_NULL);
        
        // temp1 = 1 + centrality
        LAGr_eWiseAdd(temp1, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_FP64, temp1, *centrality, GrB_NULL);
        
        // temp2 = S_matrix(i,:)'
        LAGr_extract(temp2, GrB_NULL, GrB_NULL, S_matrix, GrB_ALL, n, i, transpose_first_arg);
        
        // temp2 = (1 + centrality) ./ S_matrix(i,:)'
        LAGr_eWiseMult(temp2, GrB_NULL, GrB_NULL, GrB_DIV_FP64, temp1, temp2, GrB_NULL);
        
        // temp2 = A_matrix * (1 + centrality) ./ S_matrix(i,:)'
        LAGr_mxv(temp2, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_FP64, A_matrix, temp2, GrB_NULL);
        
        // temp1 = S_matrix(i-1,:)'
        LAGr_extract(temp1, GrB_NULL, GrB_NULL, S_matrix, GrB_ALL, n, i-1, transpose_first_arg);
        
        // Compute the entire update for the centrality scores
        // temp1 = S_matrix(i-1,:)' .* A_matrix * (1 + centrality) ./ S_matrix(i,:)'
        LAGr_eWiseMult(temp1, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_FP64, temp1, temp2, GrB_NULL);

        // Update centrality scores
        // centrality += temp1
        LAGr_assign(*centrality, GrB_NULL, GrB_PLUS_FP64, temp1, GrB_ALL, n, GrB_NULL);
    }

    // === Clean up and return =================================================
    LAGRAPH_FREE_ALL;
    return GrB_SUCCESS;
}
