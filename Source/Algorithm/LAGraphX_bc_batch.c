//------------------------------------------------------------------------------
// LAGraphX_bc_batch: Brandes' algorithm for computing betweeness centrality
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

#define LAGRAPH_FREE_WORK                   \
{                                           \
    GrB_free(&frontier);                    \
    GrB_free(&paths);                       \
    LAGraph_free(paths_dense);              \
    LAGraph_free(bc_update_dense);          \
    GrB_free(&t1);                          \
    GrB_free(&t2);                          \
    GrB_free(&desc_tsr);                    \
    GrB_free(&replace);                     \
    GrB_free(&temp);                        \
    if (S_array != NULL)                    \
    {                                       \
        for (int64_t d = 0; d < depth; d++) \
        {                                   \
            GrB_free(&(S_array[d]));        \
        }                                   \
        free (S_array);                     \
    }                                       \
}

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE_WORK ;             \
    GrB_free (centrality) ;         \
}

GrB_Info LAGraphX_bc_batch // betweeness centrality, batch algorithm
(
    GrB_Vector *centrality,    // centrality(i) is the betweeness centrality of node i
    const GrB_Matrix A_matrix, // input graph, treated as if boolean in semiring
    const GrB_Index *sources,  // source vertices from which to compute shortest paths
    int32_t num_sources        // number of source vertices (length of s)
)
{
    (*centrality) = NULL ;
    GrB_Index n; // Number of nodes in the graph

    // Array of BFS search matrices
    // S_array[i] is a matrix that stores the depth at which each vertex is 
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Matrix *S_array = NULL ;

    // Frontier matrix
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Matrix frontier = NULL ;
    
    // Paths matrix holds the number of shortest paths for each node and 
    // starting node discovered so far. Starts out sparse and becomes denser.
    GrB_Matrix paths = NULL ;
    double* paths_dense = NULL;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node. Treated as dense for efficiency.
    double* bc_update_dense = NULL;

    // Temporary workspace matrix
    GrB_Matrix temp = NULL ;

    GrB_Descriptor desc_tsr = NULL ;
    GrB_Descriptor replace = NULL ;

    GrB_Index* r = NULL;
    GrB_Index* c = NULL;
    GrB_Index* Tp = NULL;
    GrB_Index* Ti = NULL;
    double* Tx = NULL;

    GrB_Matrix t1 = NULL;
    GrB_Matrix t2 = NULL;
    double* t2_dense = NULL;

    int64_t depth = 0; // Initial BFS depth
    GxB_set (GxB_FORMAT, GxB_BY_COL) ;
    LAGr_Matrix_nrows(&n, A_matrix); // Get dimensions

    // Create a new descriptor that represents the following traits:
    //  - Tranpose the first input matrix
    //  - Replace the output 
    //  - Use the structural complement of the mask
    LAGr_Descriptor_new(&desc_tsr);
    LAGr_Descriptor_set(desc_tsr, GrB_INP0, GrB_TRAN);
    LAGr_Descriptor_set(desc_tsr, GrB_OUTP, GrB_REPLACE);
    LAGr_Descriptor_set(desc_tsr, GrB_MASK, GrB_SCMP);

    // This is also: LAGraph_desc_tocr :  A', compl mask, replace

    // Initialize paths to source vertices with ones
    // paths[s[i],i]=1 for i=[0, ..., num_sources)

    if (sources == GrB_ALL)
    {
        num_sources = n ;
    }

    LAGr_Matrix_new(&paths, GrB_INT64, n, num_sources);
    // optional: to set the matrix to CSC format
    // LAGRAPH_OK (GxB_set (paths, GxB_FORMAT, GxB_BY_COL)) ;
    if (sources == GrB_ALL)
    {
        for (GrB_Index i = 0; i < num_sources; ++i)
        {
            // paths [i,i] = 1
            LAGr_Matrix_setElement(paths, (int64_t) 1, i, i);
        }
    }
    else
    {
        for (GrB_Index i = 0; i < num_sources; ++i)
        {
            // paths [s[i],i] = 1
            LAGr_Matrix_setElement(paths, (int64_t) 1, sources[i], i);
        }
    }

    // Create frontier matrix and initialize to outgoing nodes from
    // all source nodes
    LAGr_Matrix_new(&frontier, GrB_INT64, n, num_sources);

    // AT = A'
    // frontier <!paths> = AT (:,sources)
    LAGr_extract(frontier, paths, GrB_NULL, A_matrix, GrB_ALL, n, sources, num_sources, desc_tsr);

    // Allocate memory for the array of S matrices
    S_array = (GrB_Matrix*) calloc(n, sizeof(GrB_Matrix));
    if (S_array == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    // === Breadth-first search stage ==========================================
    GrB_Index sum = 0;    // Sum of shortest paths to vertices at current depth
                          // Equal to sum(frontier). Continue BFS until new paths
                          //  are no shorter than any existing paths.
    do
    {

        // Create the current search matrix - one column for each source/BFS
        LAGr_Matrix_new(&(S_array[depth]), GrB_BOOL, n, num_sources);

        // Copy the current frontier to S
        LAGr_apply(S_array[depth], GrB_NULL, GrB_NULL, GrB_IDENTITY_BOOL, frontier, GrB_NULL);

        // Accumulate path counts: paths += frontier
        LAGr_assign(paths, GrB_NULL, GrB_PLUS_INT64, frontier, GrB_ALL, n, GrB_ALL, num_sources, GrB_NULL);
        //LAGr_eWiseAdd(paths, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_INT64, paths, frontier, GrB_NULL);

        // Update frontier: frontier<!paths>=A’ +.∗ frontier
        LAGr_mxm(frontier, paths, GrB_NULL, GxB_PLUS_TIMES_INT64, A_matrix, frontier, desc_tsr);

        // Sum up the number of BFS paths still being explored
        LAGr_Matrix_nvals(&sum, frontier);

        depth = depth + 1;

    } while (sum); // Repeat until no more shortest paths being discovered

    // === Betweenness centrality computation phase ============================

    // Create the dense update matrix and initialize it to 1
    // We will store it column-wise (col * p + row)
    const GrB_Index nnz_dense = n * num_sources;
    bc_update_dense = LAGraph_malloc(nnz_dense, sizeof(double));
    for (GrB_Index nz = 0; nz < nnz_dense; nz++)
    {
        bc_update_dense[nz] = 1.0;
    }

    // By this point, paths is (mostly) dense.
    // Create a dense version of the GraphBLAS paths matrix
    paths_dense = calloc(nnz_dense, sizeof(double));

    GrB_Index nnz_paths = nnz_dense;
    LAGr_Matrix_nvals(&nnz_paths, paths);
    GrB_Index* paths_row = LAGraph_malloc(nnz_paths, sizeof(GrB_Index));
    GrB_Index* paths_col = LAGraph_malloc(nnz_paths, sizeof(GrB_Index));
    double* paths_x = LAGraph_malloc(nnz_paths, sizeof(double));

    // NOTE: We can't use GxB_Matrix_export because paths may not be fully dense.
    LAGr_Matrix_extractTuples(paths_row, paths_col, paths_x, &nnz_paths, paths);
    for(int64_t i = 0; i < nnz_paths; i++)
    {
        // Scatter the sparse matrix values into the dense version
        GrB_Index row = paths_row[i];
        GrB_Index col = paths_col[i];
        paths_dense[col * n + row] = paths_x[i];
    }

    // Clean up our arrays
    LAGraph_free(paths_x);
    LAGraph_free(paths_row);
    LAGraph_free(paths_col);

    // Create temporary workspace matrix
    LAGr_Matrix_new(&t2, GrB_FP64, n, num_sources);

    GrB_Index nnz;
    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth - 1; i > 0; i--)
    {
        // Add contributions by successors and mask with that BFS level’s frontier

        // temp<S_array[i]> = bc_update ./ paths
        r = LAGraph_malloc(nnz_dense, sizeof(GrB_Index));
        c = LAGraph_malloc(nnz_dense, sizeof(GrB_Index));
        Tp = LAGraph_malloc(nnz_dense, sizeof(GrB_Index));
        Ti = LAGraph_malloc(nnz_dense, sizeof(GrB_Index));
        Tx = LAGraph_malloc(nnz_dense, sizeof(double));

        LAGr_Matrix_nvals(&nnz, S_array[i]);

        // Get the pattern of S_array[i]
        LAGr_Matrix_extractTuples(r, c, Tx, &nnz, S_array[i]);

        // Compute Tx = bc_update ./ paths_dense for all elements of S_array
        // Build the Tp and Ti vectors, too.
        GrB_Index nz = 0;
        GrB_Index k = 0;
        int64_t nonempty_count = 0;
        for (int64_t col = 0; col < num_sources; col++)
        {
            Tp[k++] = nz;
            bool is_nonempty_col = false;
            while (nz < nnz && col == c[nz])
            {
                // Count nonempty columns
                if (is_nonempty_col == false)
                {
                    nonempty_count += 1;
                }
                is_nonempty_col = true;

                // Compute Tx by eWiseMult of dense matrices
                GrB_Index row = Ti[nz] = r[nz];
                Tx[nz] = bc_update_dense[col * n + row] / paths_dense[col * n + row];

                nz += 1;
            }
        }
        Tp[k] = nz;

        // Create a GraphBLAS matrix t1 from Tp, Ti, Tx
        // The row/column indices are the pattern r/c from S_array[i]
        GxB_Matrix_import_CSC(&t1, GrB_FP64, n, num_sources, nz, nonempty_count, &Tp, &Ti, (void**) &Tx, GrB_NULL);

        // temp<S_array[i−1]> = (A * temp)
        LAGr_mxm(t2, S_array[i-1], GrB_NULL, GxB_PLUS_TIMES_FP64, A_matrix, t1, LAGraph_desc_ooor);

        // bc_update += t2 .* paths
        // Create a dense version of t2. Again, we can't use GxB_Matrix_export
        // because t2 may not be fully dense.
        GrB_Index nnz_t2 = nnz_dense;
        t2_dense = LAGraph_malloc(nnz_t2, sizeof(double));
        LAGr_Matrix_extractTuples(r, c, t2_dense, &nnz_t2, t2);

        for (nz = 0; nz < nnz_t2; nz++)
        {
            GrB_Index row = r[nz];
            GrB_Index col = c[nz];
            bc_update_dense[col * n + row] += t2_dense[nz] * paths_dense[col * n + row];
        }

        LAGraph_free(r);
        LAGraph_free(c);
        LAGraph_free(t2_dense);
    }

    // Initialize the centrality array with -(num_sources) to avoid counting
    // zero length paths
    double* centrality_dense = LAGraph_malloc(n, sizeof(double));
    for (int64_t i = 0; i < n; i++)
    {
        centrality_dense[i] = -num_sources;
    }

    // centrality[i] += bc_update[i,:]
    // Both are dense. We can also take care of the reduction.
    // We could also speed this up in parallel.
    for (int64_t i = 0; i < num_sources; i++)
    {
        for (int64_t j = 0; j < n; j++)
        {
            centrality_dense[j] += bc_update_dense[n * i + j];
        }
    }

    // Build the index vector.
    GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
    for (int64_t j = 0; j < n; j++)
    {
        I[j] = j;
    }

    // Import the dense vector into GraphBLAS and return it.
    GxB_Vector_import(centrality, GrB_FP64, n, n, &I, (void**) &centrality_dense, GrB_NULL);

    LAGRAPH_FREE_WORK ;
    return GrB_SUCCESS;
}
