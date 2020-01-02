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
    LAGRAPH_FREE_WORK;              \
    GrB_free (centrality);          \
}

GrB_Info LAGraphX_bc_batch3 // betweeness centrality, batch algorithm
(
    GrB_Vector *centrality,    // centrality(i): betweeness centrality of node i
    const GrB_Matrix A_matrix, // input graph
    const GrB_Matrix AT_matrix, // A'
    const GrB_Index *sources,  // source vertices for shortest paths
    int32_t num_sources        // number of source vertices (length of s)
)
{
    double tic [2];
    LAGraph_tic (tic);

    (*centrality) = NULL;
    GrB_Index n; // Number of nodes in the graph

    int nthreads ;
    GxB_get (GxB_NTHREADS, &nthreads) ;

    // Array of BFS search matrices
    // S_array[i] is a matrix that stores the depth at which each vertex is 
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Matrix *S_array = NULL;

    // Frontier matrix
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Matrix frontier = NULL;
    
    // Paths matrix holds the number of shortest paths for each node and 
    // starting node discovered so far. Starts out sparse and becomes denser.
    GrB_Matrix paths = NULL;
    double *paths_dense = NULL;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node. Treated as dense for efficiency.
    double *bc_update_dense = NULL;

    GrB_Descriptor desc_tsr = NULL;
    GrB_Descriptor replace = NULL;

    GrB_Index* Sp = NULL;
    GrB_Index* Si = NULL;
    double *Sx = NULL;

    GrB_Index* Tp = NULL;
    GrB_Index* Ti = NULL;
    double *Tx = NULL;

    GrB_Index num_rows;
    GrB_Index num_cols;
    GrB_Index nnz;
    int64_t num_nonempty;
    GrB_Type type;

    GrB_Matrix t1 = NULL;
    GrB_Matrix t2 = NULL;

    int64_t depth = 0; // Initial BFS depth

    LAGr_Matrix_nrows(&n, A_matrix); // Get dimensions

    const GrB_Index nnz_dense = n * num_sources;

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
        num_sources = n;
    }

    LAGr_Matrix_new(&paths, GrB_FP64, n, num_sources);
    GxB_set(paths, GxB_FORMAT, GxB_BY_COL);

    // make paths dense
    LAGr_assign(paths, NULL, NULL, 0, GrB_ALL, n, GrB_ALL, num_sources, NULL);

    // Force resolution of pending tuples
    GrB_Index ignore;
    GrB_Matrix_nvals(&ignore, paths);

    if (sources == GrB_ALL)
    {
        // TODO: remove this option
        for (GrB_Index i = 0; i < num_sources; ++i)
        {
            // paths [i,i] = 1
            LAGr_Matrix_setElement(paths, (double) 1, i, i);
        }
    }
    else
    {
        for (GrB_Index i = 0; i < num_sources; ++i)
        {
            // paths [s[i],i] = 1
            LAGr_Matrix_setElement(paths, (double) 1, sources[i], i);
        }
    }

    // Create frontier matrix and initialize to outgoing nodes from
    // all source nodes
    LAGr_Matrix_new(&frontier, GrB_FP64, n, num_sources);
    GxB_set(frontier, GxB_FORMAT, GxB_BY_COL);

    // AT = A'
    // frontier <!paths> = AT (:,sources)
    // TODO: use mxm, so A_matrix values are ignored.
    LAGr_extract(frontier, paths, GrB_NULL, A_matrix, GrB_ALL, n, sources, num_sources, desc_tsr);

    // Allocate memory for the array of S matrices
    S_array = (GrB_Matrix*) calloc(n, sizeof(GrB_Matrix));
    if (S_array == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL;
        return (GrB_OUT_OF_MEMORY);
    }

    printf ("starting 1st phase\n") ;
    // GxB_fprint (A_matrix, 2, stdout) ;
    // GxB_fprint (frontier, 2, stdout) ;

    //=== Breadth-first search stage ===========================================
    GrB_Index sum = 0;    // Sum of shortest paths to vertices at current depth
                          // Equal to sum(frontier). Continue BFS until new paths
                          //  are no shorter than any existing paths.

    double time_1 = LAGraph_toc (tic) ;
    printf ("Xbc setup %g sec\n", time_1) ;
    double time_2 = 0 ;
    double time_2a_all = 0 ;
    double time_2b_all = 0 ;

    do
    {
        LAGraph_tic (tic);

        // Create the current search matrix - one column for each source/BFS
        LAGr_Matrix_new(&(S_array[depth]), GrB_BOOL, n, num_sources);
        GxB_set(S_array[depth], GxB_FORMAT, GxB_BY_COL);

        // Copy the current frontier to S
        LAGr_apply(S_array[depth], GrB_NULL, GrB_NULL, GrB_IDENTITY_BOOL, frontier, GrB_NULL);

        //=== Accumulate path counts: paths += frontier ========================

        // Export paths
        GxB_Matrix_export_CSC(&paths, &type, &num_rows, &num_cols, &nnz, &num_nonempty, &Sp, &Si, (void **) &Sx, GrB_NULL);

        // Export frontier
        GxB_Matrix_export_CSC(&frontier, &type, &num_rows, &num_cols, &nnz, &num_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        // Use frontier pattern to update dense paths
#pragma omp parallel for num_threads(nthreads)
        for (int64_t col = 0; col < num_sources; col++)
        {
            for (GrB_Index p = Tp[col]; p < Tp[col+1]; p++)
            {
                GrB_Index row = Ti[p];
                Sx [col * n + row] += Tx [p];
            }
        }

        // Import frontier
        GxB_Matrix_import_CSC(&frontier, GrB_FP64, n, num_sources, nnz, num_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        // Import paths
        GxB_Matrix_import_CSC(&paths, GrB_FP64, n, num_sources, nnz_dense, n, &Sp, &Si, (void **) &Sx, GrB_NULL);

        double time_2a = LAGraph_toc (tic) ;
        time_2 += time_2a ;
        time_2a_all += time_2a ;
        // printf (" %16"PRId64"    accum: %12.6g  ", depth, time_2a) ;

        LAGraph_tic (tic);

        //=== Update frontier: frontier<!paths>=A’ +.∗ frontier ================
        // uses the "push" method (saxpy).  Note the dense complemented mask.
        LAGr_mxm(frontier, paths, GrB_NULL, GxB_PLUS_SECOND_FP64, A_matrix, frontier, desc_tsr);

        //=== Sum up the number of BFS paths still being explored ==============
        LAGr_Matrix_nvals(&sum, frontier);

        depth = depth + 1;

        double time_2b = LAGraph_toc (tic) ;
        // printf ("    mxm:   %12.6g\n", time_2b) ;
        time_2 += time_2b ;
        time_2b_all += time_2b ;

    } while (sum); // Repeat until no more shortest paths being discovered

    printf ("Xbc bfs phase: %g\n", time_2) ;
    printf ("    bfs mxm:   %g\n", time_2b_all) ;
    printf ("    bfs other  %g\n", time_2a_all) ;
    LAGraph_tic (tic);

    // GxB_fprint (paths, 3, stdout) ;

    //=== Betweenness centrality computation phase =============================

    // Create the dense update matrix and initialize it to 1
    // We will store it column-wise (col * p + row)
    bc_update_dense = LAGraph_malloc(nnz_dense, sizeof(double));
#pragma omp parallel for num_threads(nthreads)
    for (GrB_Index nz = 0; nz < nnz_dense; nz++)
    {
        bc_update_dense[nz] = 1.0;
    }

    // By this point, paths is (mostly) dense.
    // Create a dense version of the GraphBLAS paths matrix
    GxB_Matrix_export_CSC(&paths, &type, &num_rows, &num_cols, &nnz, &num_nonempty, &Sp, &Si, (void **) &paths_dense, GrB_NULL);

    // Throw away the "sparse" version of paths
    LAGraph_free(Sp);
    LAGraph_free(Si);

    // Create temporary workspace matrix
    LAGr_Matrix_new(&t2, GrB_FP64, n, num_sources);

    double time_3 = LAGraph_toc (tic) ;
    printf ("Xbc: setup for backtrack %12.6g\n", time_3) ;
    double time_4 = 0 ;
    double time_4a_all = 0 ;
    double time_4b_all = 0 ;


    printf ("starting 2nd phase\n") ;
    // GxB_fprint (AT_matrix, 2, stdout) ;

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth - 1; i > 0; i--)
    {
        // Add contributions by successors and mask with that BFS level’s frontier
        LAGraph_tic (tic);

        //=== temp<S_array[i]> = bc_update ./ paths ============================

        // Export the pattern of S_array[i]
        void *Bx ;
        GxB_Matrix_export_CSC(&(S_array[i]), &type, &num_rows, &num_cols, &nnz, &num_nonempty, &Sp, &Si, &Bx, GrB_NULL);

        // Compute Tx = bc_update ./ paths_dense for all elements of S_array
        // Build the Tp and Ti vectors, too.
        Tp = LAGraph_malloc(num_sources+1, sizeof(GrB_Index));
        Ti = LAGraph_malloc(nnz, sizeof(GrB_Index));
        Tx = LAGraph_malloc(nnz, sizeof(double));

#pragma omp parallel for num_threads(nthreads)
        for (int64_t col = 0; col < num_sources; col++)
        {
            Tp[col] = Sp[col];
            for (GrB_Index p = Sp[col]; p < Sp[col+1]; p++)
            {
                // Compute Tx by eWiseMult of dense matrices
                GrB_Index row = Ti[p] = Si[p];
                Tx [p] = bc_update_dense [col * n + row] / paths_dense[col * n + row] ;
            }
        }
        Tp[num_sources] = Sp[num_sources];

        // Restore S_array[i] by importing it
        GxB_Matrix_import_CSC(&(S_array[i]), GrB_BOOL, num_rows, num_cols, nnz, num_nonempty, &Sp, &Si, &Bx, GrB_NULL);

        // Create a GraphBLAS matrix t1 from Tp, Ti, Tx
        // The row/column indices are the pattern r/c from S_array[i]
        GxB_Matrix_import_CSC(&t1, GrB_FP64, n, num_sources, nnz, num_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        double time_4a = LAGraph_toc (tic) ;
        time_4 += time_4a ;
        time_4a_all += time_4a ;
        // printf (" %16"PRId64"    contr: %12.6g  ", i, time_4a) ;
        LAGraph_tic (tic);

        // GxB_fprint (t1, 2, stdout) ;

        //=== t2<S_array[i−1]> = (A * t1) ======================================
        // uses the "pull" method (dot) which is slow
//      LAGr_mxm(t2, S_array[i-1], GrB_NULL, GxB_PLUS_SECOND_FP64, A_matrix, t1, LAGraph_desc_ooor);
        // uses the "push" method (saxpy).  Note the sparse mask; not complemented.

        //=== t2<S_array[i−1]> = (A'' * t1) ======================================
        LAGr_mxm(t2, S_array[i-1], GrB_NULL, GxB_PLUS_SECOND_FP64, AT_matrix, t1, LAGraph_desc_toor);
        GrB_free(&t1);

        //=== bc_update += t2 .* paths =========================================
        GxB_Matrix_export_CSC(&t2, &type, &num_rows, &num_cols, &nnz, &num_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);
#pragma omp parallel for num_threads(nthreads)
        for (int64_t col = 0; col < num_sources; col++)
        {
            for (GrB_Index p = Tp[col]; p < Tp[col+1]; p++)
            {
                GrB_Index row = Ti[p];
                bc_update_dense[col * n + row] += Tx [p] * paths_dense [col * n + row] ;
            }
        }

        // Re-import t2
        GxB_Matrix_import_CSC(&t2, GrB_FP64, num_rows, num_cols, nnz, num_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        double time_4b = LAGraph_toc (tic) ;
        time_4 += time_4b ;
        time_4b_all += time_4b ;
        // printf ("    mxm:   %12.6g\n", time_4b) ;

    }

    printf ("Xbx 2nd phase: %g\n", time_4) ;
    printf ("    2nd mxm  : %g\n", time_4b_all) ;
    printf ("    2nd other: %g\n", time_4a_all) ;
    LAGraph_tic (tic);

    //=== Initialize the centrality array with -(num_sources) to avoid counting
    //    zero length paths ====================================================
    double *centrality_dense = LAGraph_malloc(n, sizeof(double));
#pragma omp parallel for num_threads(nthreads)
    for (GrB_Index i = 0; i < n; i++)
    {
        centrality_dense[i] = -num_sources;
    }

    //=== centrality[i] += bc_update[i,:] ======================================
    // Both are dense. We can also take care of the reduction.
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (GrB_Index j = 0; j < n; j++)
    {
        for (int64_t i = 0; i < num_sources; i++)
        {
            centrality_dense[j] += bc_update_dense[n * i + j];
        }
    }

    // Build the index vector.
    GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
#pragma omp parallel for num_threads(nthreads)
    for (GrB_Index j = 0; j < n; j++)
    {
        I[j] = j;
    }

    // Import the dense vector into GraphBLAS and return it.
    GxB_Vector_import(centrality, GrB_FP64, n, n, &I, (void **) &centrality_dense, GrB_NULL);

    LAGRAPH_FREE_WORK;
    double time_5 = LAGraph_toc (tic) ;
    printf ("Xbc wrapup:    %g\n", time_5) ;
    printf ("Xbc total:     %g\n", time_1 + time_2 + time_3 + time_4 + time_5) ;

    return GrB_SUCCESS;
}
