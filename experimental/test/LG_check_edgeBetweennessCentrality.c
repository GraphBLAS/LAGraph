//------------------------------------------------------------------------------
// LG_check_edgeBetweennessCentrality: reference implementation for edge 
// betweenness centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Casey Pei, Texas A&M University

//------------------------------------------------------------------------------

#define LG_FREE_WORK                                \
{                                                   \
    LAGraph_Free ((void **) &queue, NULL) ;         \
    LAGraph_Free ((void **) &depth, NULL) ;         \
    LAGraph_Free ((void **) &bc_vertex_flow, NULL) ;\
    LAGraph_Free ((void **) &S, NULL) ;             \
    LAGraph_Free ((void **) &paths, NULL) ;         \
    LAGraph_Free ((void **) &Pj, NULL) ;            \
    LAGraph_Free ((void **) &Ptail, NULL) ;         \
    if(AT != G->AT)                                 \
    {                                               \
        GrB_free (&AT) ;                            \
    }                                               \
}

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    LAGraph_Free ((void **) &Ap, NULL) ;            \
    LAGraph_Free ((void **) &Aj, NULL) ;            \
    LAGraph_Free ((void **) &Ax, NULL) ;            \
    LAGraph_Free ((void **) &ATp, NULL) ;           \
    LAGraph_Free ((void **) &ATj, NULL) ;           \
    LAGraph_Free ((void **) &ATx, NULL) ;           \
    LAGraph_Free ((void **) &result, NULL) ;        \
}

#include "LG_internal.h"
#include <LAGraphX.h>
#include "LG_Xtest.h"


//------------------------------------------------------------------------------
// test the results from a Edge Betweenness Centrality
//------------------------------------------------------------------------------

int LG_check_edgeBetweennessCentrality
(
    // output
    GrB_Matrix *C,      // centrality matrix
    // input
    LAGraph_Graph G,
    GrB_Vector sources,         // source vertices to compute shortest paths (if NULL or empty, use all vertices)
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // initialize workspace variables
    //--------------------------------------------------------------------------

    double tt = LAGraph_WallClockTime ( ) ;
    GrB_Info info ;

    // Array storing shortest path distances from source to each vertex
    int64_t *depth = NULL ;

    // Array storing dependency scores during accumulation phase
    double *bc_vertex_flow = NULL ;

    // Stack used for backtracking phase in dependency accumulation
    int64_t *S = NULL ;

    // Queue used for BFS traversal
    int64_t *queue = NULL ;

    // Predecessor list components:
    // Pj: array of predecessor vertices
    // Ptail: end indices for each vertex's predecessor list
    // Phead: start indices for each vertex's predecessor list
    GrB_Index *Pj = NULL ;
    GrB_Index *Ptail = NULL ;
    GrB_Index *Phead = NULL ;

    // Array storing number of shortest paths to each vertex
    double *paths = NULL ;

    // Temporary array for centrality results
    double *result = NULL;

    GrB_Vector internal_sources = NULL;
    bool created_sources = false;

    GrB_Matrix AT  = NULL ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Index *Ap = NULL, *Aj = NULL, *neighbors = NULL, *ATp = NULL, *ATj = NULL ;
    void *Ax = NULL, *ATx = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, nvals, ATp_size, ATj_size, ATx_size ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
    #if defined ( COVERAGE )
    bool print_timings = true ;
    #else
    bool print_timings = (n >= 2000) ;
    #endif

    LG_TRY (LAGraph_DeleteSelfEdges (G, msg)) ;

    GrB_Matrix A = G->A ;

    LG_TRY (LAGraph_Cached_AT (G, msg)) ;

    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
         G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        // AT = A;
        // GrB_Matrix_new (&AT, GrB_FP64, n, n) ;
        GrB_Matrix_dup (&AT, A) ;
    }
    else
    {
        // A and A' differ
         AT = G->AT ;
         LG_ASSERT_MSG (AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required") ;
    }

    
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    
    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    LG_TRY(LAGraph_Malloc((void **)&depth, n, sizeof(int64_t), msg));

    LG_TRY(LAGraph_Calloc((void **)&bc_vertex_flow, n, sizeof(double), msg));

    LG_TRY(LAGraph_Malloc((void **)&S, n, sizeof(int64_t), msg));

    LG_TRY(LAGraph_Malloc((void **)&queue, n, sizeof(int64_t), msg));

    //--------------------------------------------------------------------------
    // bfs on the A
    //--------------------------------------------------------------------------

    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("LG_check_bfs init  time: %g sec\n", tt) ;
        tt = LAGraph_WallClockTime ( ) ;
    }

    // Initialize centrality matrix result to 0
    // 1. result [(v, w)] ← 0, ∀(v, w) ∈ E
    // A temporary result centrality matrix initialized to 0 for all vertice,
    // -- further changes would need to be made to make it a dictionary of edges.
    GrB_Index result_size = n * n ;
    LG_TRY(LAGraph_Calloc((void **)&result, result_size, sizeof(double), msg));

    // result (v,w) is held in result (INDEX(v,w)):
    #define INDEX(i,j) ((i)*n+(j))

    //--------------------------------------------------------------------------
    // unpack the A matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    bool iso, AT_iso ; 
    GRB_TRY (GxB_Matrix_unpack_CSR (A,
        &Ap, &Aj, &Ax, &Ap_size, &Aj_size, &Ax_size, &iso, NULL, NULL)) ;

    GRB_TRY (GxB_Matrix_unpack_CSR (AT,
        &ATp, &ATj, &ATx, &ATp_size, &ATj_size, &ATx_size, &AT_iso, NULL, NULL)) ;
    #endif

    Phead = ATp ;

    //--------------------------------------------------------------------------

    LG_TRY(LAGraph_Malloc((void **)&Pj, nvals, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&Ptail, n, sizeof(GrB_Index), msg)); // might need to be + 1

    LAGraph_Calloc ((void **) &paths, n, sizeof (double), msg) ;

    if (sources == NULL)
    {
        // Create a vector with all nodes as sources
        GRB_TRY (GrB_Vector_new (&internal_sources, GrB_INT64, n)) ;

        // internal_sources (0:n-1) = 0
        GRB_TRY (GrB_assign (internal_sources, NULL, NULL, 0, GrB_ALL, n, NULL)) ;

        // internal_sources (0:n-1) = 0:n-1
        GRB_TRY (GrB_apply (internal_sources, NULL, NULL, GrB_ROWINDEX_INT64,
            internal_sources, 0, NULL)) ;
        
        // Use this vector instead
        sources = internal_sources;
        created_sources = true;
    }
    
    // Extract number of source nodes
    GrB_Index nvals_sources;
    GRB_TRY (GrB_Vector_nvals (&nvals_sources, sources));
    
    if (nvals_sources == 0)
    {
        // If sources vector is empty, return an empty centrality matrix
        // (Create an empty C matrix or set it to zeros as needed)
        GRB_TRY (GrB_Matrix_new(C, GrB_FP64, n, n));
        
        // Clean up resources
        if (created_sources) GRB_TRY (GrB_free(&internal_sources));
        LG_FREE_ALL ;
        return (GrB_SUCCESS);
    }

    // =========================================================================
    // === Main computation loop ==============================================
    // =========================================================================

    GrB_Index s;

    // Process each source vertex
    for (GrB_Index i = 0; i < nvals_sources; i++) {
        GRB_TRY (GrB_Vector_extractElement(&s, sources, i)) ;

        // check for invalid indices
        LG_ASSERT (s < n, GrB_INVALID_VALUE) ;

        size_t sp = 0;  // stack pointer
        memcpy(Ptail, ATp, n * sizeof(GrB_Index));

        // Initialize path counts
        for (int64_t i = 0; i < n; i++) {
            paths[i] = 0;
        }
        paths[s] = 1;

        // Initialize distances
        for (size_t t = 0; t < n; t++) {
            depth[t] = -1;
        }
        depth[s] = 0;

        //----------------------------------------------------------------------
        // BFS phase to compute shortest paths
        //----------------------------------------------------------------------

        int64_t qh = 0, qt = 0;  // queue head and tail
        queue[qt++] = s;         // enqueue source

        while (qh < qt) {
            int64_t v = queue[qh++];
            S[sp++] = v;

            // Process neighbors of current vertex
            for (int64_t p = Ap[v]; p < Ap[v+1]; p++) {
                int64_t w = Aj[p];
                
                // Handle unvisited vertices
                if (depth[w] < 0) {
                    queue[qt++] = w;
                    depth[w] = depth[v] + 1;
                }

                // Update path counts for vertices at next level
                if (depth[w] == depth[v] + 1) {
                    paths[w] += paths[v];

#if 1
                    LG_ASSERT (Ptail [w] < Phead [w+1], GrB_INVALID_VALUE) ;
                    LG_ASSERT (Ptail [w] >= Phead [w], GrB_INVALID_VALUE) ;
#else
                    if (Ptail [w] >= Phead [w+1] || Ptail [w] < Phead [w])
                    {
                        printf ("Ack! w=%ld Ptail [w]=%ld, Phead [w]=%ld Phead[w+1]=%ld\n", 
                            w, Ptail [w], Phead [w], Phead [w+1]) ;
                        fflush (stdout) ; abort ( ) ;
                    }
#endif

                    Pj[Ptail[w]++] = v;

                }
            }   
        }

        //----------------------------------------------------------------------
        // Dependency accumulation phase
        //----------------------------------------------------------------------

        // Initialize dependency scores
        for (size_t v = 0; v < n; v++) {
            bc_vertex_flow[v] = 0;
        }

        // Process vertices in reverse order of discovery
        while (sp > 0) {
            int64_t w = S[--sp];

            // Update dependencies through predecessors
            for (int64_t p = Phead[w]; p < Ptail[w]; p++) {
                int64_t v = Pj[p];

                // Compute and accumulate dependency
                double centrality = paths[v] * ((bc_vertex_flow[w] + 1) / paths[w]);
                bc_vertex_flow[v] += centrality;

                if (G->kind == LAGraph_ADJACENCY_UNDIRECTED) {
                    result[INDEX(v,w)] += centrality / 2;
                    result[INDEX(w,v)] += centrality / 2;
                }
                else {
                    result[INDEX(v,w)] += centrality;
                }
            }
        }
    }

    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("LG_check_edgeBetweenessCentrality time: %g sec\n", tt) ;
        tt = LAGraph_WallClockTime ( ) ;
    }

    //--------------------------------------------------------------------------
    // repack the A matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_Matrix_pack_CSR (A,
        &Ap, &Aj, &Ax, Ap_size, Aj_size, Ax_size, iso, false, NULL)) ;
    GRB_TRY (GxB_Matrix_pack_CSR (AT,
        &ATp, &ATj, &ATx, ATp_size, ATj_size, ATx_size, AT_iso, false, NULL)) ;
    #endif

#if 0
GrB_Info GxB_Matrix_pack_FullR  // pack a full matrix, held by row
(
    GrB_Matrix A,       // matrix to create (type, nrows, ncols unchanged)
    void **Ax,          // values, Ax_size >= nrows*ncols * (type size)
                        // or Ax_size >= (type size), if iso is true
    GrB_Index Ax_size,  // size of Ax in bytes
    bool iso,           // if true, A is iso
    const GrB_Descriptor desc
) ;
#endif

    GrB_Matrix C_temp;
    LG_TRY (GrB_Matrix_new(&C_temp, GrB_FP64, n, n)) ;
    LG_TRY (GxB_Matrix_pack_FullR(C_temp, (void **) &result, result_size * sizeof(double), false, NULL) ) ;

    LG_TRY (GrB_assign(C_temp, A, NULL, C_temp, GrB_ALL, n, GrB_ALL, n, GrB_DESC_RS)) ;

    *C = C_temp;

    if (created_sources) {
        GRB_TRY (GrB_free(&internal_sources));
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------
    LG_FREE_WORK ;

    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("LG_check_edgeBetweennessCentrality check time: %g sec\n", tt) ;
    }
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
