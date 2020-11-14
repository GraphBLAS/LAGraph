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

#define DO_PULL 0

// If DO_PULL is #defined, the algorithm does each GrB_mxm twice: once with the
// "pull" (dot product method internally in GraphBLAS) and ones with the "push"
// (the saxpy method in GraphBLAS).  Then it pretends to have a perfect
// heuristic by taking the min of both times to compute the "pushpull" time.
// This is of course unrealistic, but it's a lower bound on any heuristc that
// tries to select the correct method at each step.

//------------------------------------------------------------------------------

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_WORK                       \
{                                               \
    GrB_free(&frontier);                        \
    GrB_free(&paths);                           \
    LAGraph_free(paths_dense);                  \
    LAGraph_free(bc_update_dense);              \
    GrB_free(&t1);                              \
    GrB_free(&t2);                              \
    GrB_free (&pull_descriptor) ;               \
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
    LAGRAPH_FREE_WORK;              \
    GrB_free (centrality);          \
}

// TODO add LAGraph_PLUS_SECOND_FP* to LAGraph.h.
#if 0
// select FP64
#define REAL_t                  double
#define LAGr_REAL_TYPE          GrB_FP64
#define LAGr_PLUS_SECOND_REAL   GxB_PLUS_SECOND_FP64
#else
// select FP32
#define REAL_t                  float
#define LAGr_REAL_TYPE          GrB_FP32
#define LAGr_PLUS_SECOND_REAL   GxB_PLUS_SECOND_FP32
#endif

GrB_Info LAGraphX_bc_batch3 // betweeness centrality, batch algorithm
(
    GrB_Vector *centrality,    // centrality(i): betweeness centrality of node i
    const GrB_Matrix A_matrix, // input graph
    const GrB_Matrix AT_matrix, // A'
    const GrB_Index *sources,  // source vertices for shortest paths
    int32_t num_sources,                 // number of source vertices (length of s)
    double timing [3]
)
{
    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    return (GrB_NO_VALUE) ;
    #else

    GrB_Info info ;
    GrB_Descriptor pull_descriptor = NULL ;

    // Frontier matrix
    // Stores # of shortest paths to vertices at current BFS depth
    GrB_Matrix frontier = NULL;

    // Array of BFS search matrices
    // S_array[i] is a matrix that stores the depth at which each vertex is
    // first seen thus far in each BFS at the current depth i. Each column
    // corresponds to a BFS traversal starting from a source node.
    GrB_Matrix *S_array = NULL;

    // Paths matrix holds the number of shortest paths for each node and
    // starting node discovered so far. Starts out sparse and becomes denser.
    GrB_Matrix paths = NULL;
    REAL_t *paths_dense = NULL;

    // Update matrix for betweenness centrality, values for each node for
    // each starting node. Treated as dense for efficiency.
    REAL_t *bc_update_dense = NULL;

    GrB_Matrix t1 = NULL;
    GrB_Matrix t2 = NULL;

    GrB_Index n; // Number of nodes in the graph

    (*centrality) = NULL;

    //--------------------------------------------------------------------------

    double tic [2];
    LAGraph_tic (tic);

    GxB_Format_Value a_fmt, at_fmt ;
    LAGRAPH_OK (GxB_get (A_matrix,  GxB_FORMAT, &a_fmt )) ;
    LAGRAPH_OK (GxB_get (AT_matrix, GxB_FORMAT, &at_fmt)) ;
    if (a_fmt != GxB_BY_ROW || at_fmt != GxB_BY_ROW)
    {
        LAGRAPH_ERROR ("A and AT must be stored by row", GrB_INVALID_VALUE) ;
    }

    int nthreads ;
    GxB_get (GxB_NTHREADS, &nthreads) ;

    GrB_Index* Sp = NULL;
    GrB_Index* Si = NULL;
    REAL_t *Sx = NULL;

    GrB_Index* Tp = NULL;
    GrB_Index* Ti = NULL;
    REAL_t *Tx = NULL;

    GrB_Index num_rows, num_cols, nnz, anz ;
    GrB_Type type ;

    LAGr_Matrix_nrows (&n, A_matrix) ;          // # of nodes
    LAGr_Matrix_nvals (&anz, A_matrix) ;        // # of edges
    double d = ((double) anz) / ((double) n) ;  // average degree

    // descriptor for "pull" method: LAGraph_desc_oocr + dot
    GrB_Descriptor_new (&pull_descriptor) ;
    GrB_Descriptor_set (pull_descriptor, GrB_MASK, GrB_SCMP) ;
    GrB_Descriptor_set (pull_descriptor, GrB_OUTP, GrB_REPLACE) ;
    GrB_Descriptor_set (pull_descriptor, GxB_AxB_METHOD, GxB_AxB_DOT) ;

    // Initialize paths to source vertices with ones
    // paths[s[i],i]=1 for i=[0, ..., num_sources)

    if (sources == GrB_ALL)
    {
        num_sources = n;        // TODO delete this option
    }

    const GrB_Index nnz_dense = n * num_sources ;
    double ns = num_sources ;

    LAGr_Matrix_new(&paths, LAGr_REAL_TYPE, n, num_sources);
    GxB_set(paths, GxB_FORMAT, GxB_BY_COL);

    // make paths dense
    LAGr_assign (paths, NULL, NULL, 0, GrB_ALL, n, GrB_ALL, num_sources, NULL) ;

    // Force resolution of pending tuples
    GrB_Index ignore;
    GrB_Matrix_nvals(&ignore, paths);

    if (sources == GrB_ALL)
    {
        // TODO: remove this option
        for (GrB_Index i = 0; i < num_sources; ++i)
        {
            // paths [i,i] = 1
            LAGr_Matrix_setElement(paths, (REAL_t) 1, i, i);
        }
    }
    else
    {
        for (GrB_Index i = 0; i < num_sources; ++i)
        {
            // paths [s[i],i] = 1
            LAGr_Matrix_setElement(paths, (REAL_t) 1, sources[i], i);
        }
    }

    // Create frontier matrix and initialize to outgoing nodes from
    // all source nodes
    LAGr_Matrix_new(&frontier, LAGr_REAL_TYPE, n, num_sources);
    GxB_set(frontier, GxB_FORMAT, GxB_BY_COL);

    // AT = A'
    // frontier <!paths> = AT (:,sources)
    // TODO: use mxm, so A_matrix values are ignored.
    LAGr_extract(frontier, paths, GrB_NULL, A_matrix, GrB_ALL, n, sources,
        num_sources, LAGraph_desc_tocr);

    // Allocate memory for the array of S matrices
    S_array = (GrB_Matrix*) LAGraph_calloc (n, sizeof(GrB_Matrix));
    if (S_array == NULL)
    {
        // out of memory
        LAGRAPH_FREE_ALL;
        return (GrB_OUT_OF_MEMORY);
    }

    //--------------------------------------------------------------------------
    // Breadth-first search stage
    //--------------------------------------------------------------------------

    GrB_Index frontier_size = 0 ;   // size of current frontier
    LAGr_Matrix_nvals (&frontier_size, frontier) ;
    GrB_Index seen = 0 ;            // total # of nodes seen * (# sources)

    double time_1 = LAGraph_toc (tic) ;
    // printf ("   init:      %g\n", time_1) ;
    double phase1_other_time = 0 ;
    double phase1_allpush_time = 0 ;
    double phase1_allpull_time = 0 ;
    double phase1_pushpull_time = 0 ;

    int nth = LAGRAPH_MIN (nthreads, num_sources) ;

    int64_t depth = 0; // Initial BFS depth

    do
    {
        LAGraph_tic (tic);
        // printf ("depth: %g\n", (double) depth) ;

        // Create the current search matrix - one column for each source/BFS
        LAGr_Matrix_new (&(S_array[depth]), GrB_BOOL, n, num_sources) ;
        GxB_set (S_array[depth], GxB_FORMAT, GxB_BY_COL) ; 

        // Copy the current frontier to S
        LAGr_apply (S_array[depth], GrB_NULL, GrB_NULL, GrB_IDENTITY_BOOL,
            frontier, GrB_NULL) ;

        //=== Accumulate path counts: paths += frontier ========================

        // Export paths
        int64_t paths_nonempty ;
        GxB_Matrix_export_CSC(&paths, &type, &num_rows, &num_cols, &nnz,
            &paths_nonempty, &Sp, &Si, (void **) &Sx, GrB_NULL);

        // Export frontier
        int64_t frontier_nonempty ;
        GxB_Matrix_export_CSC(&frontier, &type, &num_rows, &num_cols, &nnz,
            &frontier_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        // Use frontier pattern to update dense paths
        #pragma omp parallel for num_threads(nth)
        for (int64_t col = 0; col < num_sources; col++)
        {
            for (GrB_Index p = Tp[col]; p < Tp[col+1]; p++)
            {
                GrB_Index row = Ti[p];
                Sx [col * n + row] += Tx [p];
            }
        }

        // Import frontier
        GxB_Matrix_import_CSC(&frontier, LAGr_REAL_TYPE, n, num_sources, nnz,
            frontier_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        // Import paths
        GxB_Matrix_import_CSC(&paths, LAGr_REAL_TYPE, n, num_sources,
            nnz_dense, paths_nonempty, &Sp, &Si, (void **) &Sx, GrB_NULL);

        phase1_other_time += LAGraph_toc (tic) ;

        //=== Update frontier: frontier<!paths>=A’ +.∗ frontier ================

        seen += frontier_size ;         // # nonzeros in paths array

        /*
        double u = (n * ns - seen) ;    // # zeros in paths array
        double f = frontier_size / ns ;
        double push_work_estimate = d * frontier_size ;
        double pull_work_estimate = u * fmin (d + f, d * log2 (f)) ;
        printf ("\n1: d: %g f: %g u: %g ns: %g "
            "push: %g pull: %g pull/push %g\n", d, f, u, ns,
            push_work_estimate, pull_work_estimate,
            pull_work_estimate / push_work_estimate) ;
        */

            double pull_time = INFINITY ;
#if DO_PULL
            GrB_Matrix frontier2  = NULL ;
            GrB_Matrix_dup (&frontier2, frontier) ;
            // uses the "pull" method (dot), because AT_matrix is stored by
            // row, and frontier is stored by column.
            LAGraph_tic (tic);
            LAGr_mxm(frontier2, paths, GrB_NULL, LAGr_PLUS_SECOND_REAL,
                AT_matrix, frontier2, pull_descriptor) ;
            pull_time = LAGraph_toc (tic) ;
            // printf ("1: pull_time: %g sec\n", pull_time) ;
            GrB_free (&frontier2) ;
#endif
            phase1_allpull_time += pull_time ;

            // uses the "push" method (saxpy)
            LAGraph_tic (tic);
            LAGr_mxm(frontier, paths, GrB_NULL, LAGr_PLUS_SECOND_REAL,
                A_matrix, frontier, LAGraph_desc_tocr);

            double push_time = LAGraph_toc (tic) ;
            // printf ("1: push_time: %g sec,  pull/push %g\n",
            //     push_time, pull_time/push_time) ;
            phase1_allpush_time += push_time ;

            // assume a perfect pushpull heuristic
            double pushpull_time = fmin (pull_time, push_time) ;
            phase1_pushpull_time += pushpull_time ;

        //=== Find the new frontier size =======================================
        LAGraph_tic (tic);
        LAGr_Matrix_nvals (&frontier_size, frontier) ;
        depth = depth + 1;
        phase1_other_time += LAGraph_toc (tic) ;

    } while (frontier_size > 0) ; // Repeat until the frontier is empty

    // printf ("    1st mxm allpush:  %g\n", phase1_allpush_time) ;
#if DO_PULL
    printf ("    1st mxm allpull:  %g\n", phase1_allpull_time) ;
    printf ("    1st mxm pushpull: %g\n", phase1_pushpull_time) ;
#endif
    // printf ("    1st other:        %g\n", phase1_other_time) ;

    LAGraph_tic (tic);

    //--------------------------------------------------------------------------
    // Betweenness centrality computation phase
    //--------------------------------------------------------------------------

    // Create the dense update matrix and initialize it to 1
    // We will store it column-wise (col * p + row)
    bc_update_dense = LAGraph_malloc(nnz_dense, sizeof(REAL_t));

    #pragma omp parallel for num_threads(nthreads)
    for (GrB_Index nz = 0; nz < nnz_dense; nz++)
    {
        bc_update_dense[nz] = 1.0;
    }

    // By this point, paths is (mostly) dense.
    // Create a dense version of the GraphBLAS paths matrix
    int64_t paths_nonempty ;
    GxB_Matrix_export_CSC(&paths, &type, &num_rows, &num_cols, &nnz,
        &paths_nonempty, &Sp, &Si, (void **) &paths_dense, GrB_NULL);

    // Throw away the "sparse" version of paths
    LAGraph_free(Sp);
    LAGraph_free(Si);

    // Create temporary workspace matrix
    LAGr_Matrix_new(&t2, LAGr_REAL_TYPE, n, num_sources);
    GxB_set(t2, GxB_FORMAT, GxB_BY_COL);

    double time_3 = LAGraph_toc (tic) ;
    double phase2_other_time = 0 ;
    double phase2_allpush_time = 0 ;
    double phase2_allpull_time = 0 ;
    double phase2_pushpull_time = 0 ;

    // Backtrack through the BFS and compute centrality updates for each vertex
    for (int64_t i = depth - 1; i > 0; i--)
    {
        // Add contributions by successors and mask with that BFS level's
        // frontier
        LAGraph_tic (tic);
        // printf ("back: %g\n", (double) i) ;

        /*
        GrB_Index prior_size ;
        GrB_Matrix_nvals (&prior_size,    S_array [i-1]) ;
        GrB_Matrix_nvals (&frontier_size, S_array [i]) ;
        double u = prior_size ;         // # entries in the mask
        double f = frontier_size / ns ;
        double push_work_estimate = d * frontier_size ;
        double pull_work_estimate = u * fmin (d + f, d * log2 (f)) ;
        printf ("\n2: d: %g f: %g u: %g ns: %g "
            "push: %g pull: %g pull/push %g\n", d, f, u, ns,
            push_work_estimate, pull_work_estimate,
            pull_work_estimate / push_work_estimate) ;
        */

        //=== temp<S_array[i]> = bc_update ./ paths ============================

        // Export the pattern of S_array[i]
        void *Bx ;
        int64_t S_nonempty ;
        GxB_Matrix_export_CSC(&(S_array[i]), &type, &num_rows, &num_cols, &nnz,
            &S_nonempty, &Sp, &Si, &Bx, GrB_NULL);

        // Compute Tx = bc_update ./ paths_dense for all elements of S_array
        // Build the Tp and Ti vectors, too.
        Tp = LAGraph_malloc(num_sources+1, sizeof(GrB_Index));
        Ti = LAGraph_malloc(nnz, sizeof(GrB_Index));
        Tx = LAGraph_malloc(nnz, sizeof(REAL_t));

        #pragma omp parallel for num_threads(nthreads)
        for (int64_t col = 0; col < num_sources; col++)
        {
            Tp[col] = Sp[col];
            for (GrB_Index p = Sp[col]; p < Sp[col+1]; p++)
            {
                // Compute Tx by eWiseMult of dense matrices
                GrB_Index row = Ti[p] = Si[p];
                Tx [p] = bc_update_dense [col * n + row]
                           / paths_dense [col * n + row] ;
            }
        }
        Tp[num_sources] = Sp[num_sources];

        // Restore S_array[i] by importing it
        GxB_Matrix_import_CSC(&(S_array[i]), GrB_BOOL, num_rows, num_cols,
            nnz, S_nonempty, &Sp, &Si, &Bx, GrB_NULL);

        // Create a GraphBLAS matrix t1 from Tp, Ti, Tx
        // The row/column indices are the pattern r/c from S_array[i]
        GxB_Matrix_import_CSC(&t1, LAGr_REAL_TYPE, n, num_sources, nnz,
            S_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        phase2_other_time += LAGraph_toc (tic) ;

        //=== t2<S_array[i−1]> = (A * t1) ======================================

            double pull_time = INFINITY ;
#if DO_PULL
            // uses the "pull" method (dot)
            LAGraph_tic (tic);
            GrB_free (&t2) ;
            LAGr_Matrix_new(&t2, LAGr_REAL_TYPE, n, num_sources);
            GxB_set(t2, GxB_FORMAT, GxB_BY_COL);
            LAGr_mxm(t2, S_array[i-1], GrB_NULL, LAGr_PLUS_SECOND_REAL,
                A_matrix, t1, LAGraph_desc_ooor);
            pull_time = LAGraph_toc (tic) ;
            printf ("2: pull_time: %g sec\n", pull_time) ;
#endif
            phase2_allpull_time += pull_time ;

            // uses the "push" method (saxpy)
            LAGraph_tic (tic);
            GrB_free (&t2) ;
            LAGr_Matrix_new(&t2, LAGr_REAL_TYPE, n, num_sources);
            GxB_set(t2, GxB_FORMAT, GxB_BY_COL);
            LAGr_mxm(t2, S_array[i-1], GrB_NULL, LAGr_PLUS_SECOND_REAL,
                AT_matrix, t1, LAGraph_desc_toor);

            double push_time = LAGraph_toc (tic) ;
            // printf ("2: push_time: %g sec,  pull/push %g\n", push_time,
                // pull_time/push_time) ;
            phase2_allpush_time += push_time ;

            // assume a perfect pushpull heuristic
            double pushpull_time = fmin (pull_time, push_time) ;
            phase2_pushpull_time += pushpull_time ;

        LAGraph_tic (tic);
        GrB_free(&t1);

        //=== bc_update += t2 .* paths =========================================
        int64_t t2_nonempty ;
        GxB_Matrix_export_CSC(&t2, &type, &num_rows, &num_cols, &nnz,
            &t2_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        #pragma omp parallel for num_threads(nth)
        for (int64_t col = 0; col < num_sources; col++)
        {
            for (GrB_Index p = Tp[col]; p < Tp[col+1]; p++)
            {
                GrB_Index row = Ti[p];
                bc_update_dense [col * n + row] += 
                    Tx [p] * paths_dense [col * n + row] ;
            }
        }

        // Re-import t2
        GxB_Matrix_import_CSC(&t2, LAGr_REAL_TYPE, num_rows, num_cols, nnz,
            t2_nonempty, &Tp, &Ti, (void **) &Tx, GrB_NULL);

        phase2_other_time += LAGraph_toc (tic) ;
    }

    // printf ("    2nd mxm allpush:  %g\n", phase2_allpush_time) ;
#if DO_PULL
    printf ("    2nd mxm allpull:  %g\n", phase2_allpull_time) ;
    printf ("    2nd mxm pushpull: %g\n", phase2_pushpull_time) ;
#endif
    // printf ("    2nd other:        %g\n", phase2_other_time + time_3) ;

    LAGraph_tic (tic);

    //--------------------------------------------------------------------------
    // finalize centrality scores
    //--------------------------------------------------------------------------

    //=== Initialize the centrality array with -(num_sources) to avoid counting
    //    zero length paths ====================================================
    REAL_t *centrality_dense = LAGraph_malloc(n, sizeof(REAL_t));

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

    #if GxB_IMPLEMENTATION >= GxB_VERSION (4,0,0)
    GxB_Vector_import_Full (centrality, LAGr_REAL_TYPE, n,
        (void **) &centrality_dense, GrB_NULL) ;
    #else
    // Build the index vector.
    GrB_Index* I = LAGraph_malloc(n, sizeof(GrB_Index));
    
    #pragma omp parallel for num_threads(nthreads)
    for (GrB_Index j = 0; j < n; j++)
    {
        I[j] = j;
    }

    // Import the dense vector into GraphBLAS and return it.
    GxB_Vector_import(centrality, LAGr_REAL_TYPE, n, n, &I,
        (void **) &centrality_dense, GrB_NULL);
    #endif

    LAGRAPH_FREE_WORK;
    double time_5 = LAGraph_toc (tic) ;
    // printf ("   wrapup:    %g\n", time_5) ;

    timing [0] = time_1
              + (phase1_pushpull_time + phase1_other_time)
              + time_3 + (phase2_pushpull_time + phase2_other_time) + time_5 ;

    timing [1] = time_1
                + (phase1_allpush_time  + phase1_other_time)
                + time_3 + (phase2_allpush_time  + phase2_other_time) + time_5 ;

    timing [2] = time_1
                + (phase1_allpull_time  + phase1_other_time)
                + time_3 + (phase2_allpull_time  + phase2_other_time) + time_5 ;

#if DO_PULL
    printf ("Xbc total (pushpull):    %g\n", timing [0]) ;
#endif
    // printf ("Xbc total (allpush):     %g\n", timing [1]) ;
#if DO_PULL
    printf ("Xbc total (allpull):     %g\n", timing [2]) ;
#endif

    return GrB_SUCCESS;
    #endif
}

