//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_edgeBetweennessCentrality: test for Edge
// Betweenness Centrality
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
    free(queue) ;                                   \
    free(sigma) ;                                   \
    free(d) ;                                       \
    free(delta) ;                                   \
    free(S) ;                                       \
    free(P) ;                                       \
}

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    LAGraph_Free ((void **) &Ap, NULL) ;            \
    LAGraph_Free ((void **) &Aj, NULL) ;            \
    LAGraph_Free ((void **) &Ax, NULL) ;            \
}

#include "LG_internal.h"
#include "LG_test.h"

//------------------------------------------------------------------------------
// test the results from a Edge Betweenness Centrality
//------------------------------------------------------------------------------

int test_edgeBetweenessCentrality
(
    // output
    GrB_Matrix *C,      // centrality matrix
    // input
    LAGraph_Graph G,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    double tt = LAGraph_WallClockTime ( ) ;

    GrB_Index *Ap = NULL, *Aj = NULL, *neighbors = NULL ;
    void *Ax = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, ncols ;
    int64_t *queue = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    bool print_timings = (n >= 2000) ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    GrB_Matrix A = G->A ;
    GrB_Matrix AT ;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        AT = A ;
    }
    else
    {
        // A and A' differ
        AT = G->AT ;
        LG_ASSERT_MSG (AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required") ;
    }

    LG_CLEAR_MSG ;

    //--------------------------------------------------------------------------
    // initialize / allocate workspace
    //--------------------------------------------------------------------------

    GrB_Info info;

    // Holds the distances (depth levels) from the source vertex.
    int64_t *d = malloc(n * sizeof(int64_t)) ;

    // Stores dependency scores for each vertex.
    int64_t *delta = malloc(n * sizeof(int64_t)) ;

    int64_t *S = malloc(n * sizeof(int64_t)) ;

    int64_t *queue = malloc(n * sizeof(int64_t)) ;

    //--------------------------------------------------------------------------
    // unpack the A matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    bool iso ;
    GRB_TRY (GxB_Matrix_unpack_CSR (G->A,
        &Ap, &Aj, &Ax, &Ap_size, &Aj_size, &Ax_size, &iso, NULL, NULL)) ;
    #endif

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
    // TODO make this a copy of A except with 1 = 0
    // A temporary result centrality matrix initialized to 0 for all vertice,
    // -- further changes would need to be made to make it a dictionary of edges.
    int64_t *result = calloc(n * n, sizeof(int64_t));
    result_p = malloc (n * sizeof (int64_t)) ;
    result_j = malloc (Aj_size * sizeof (int64_t)) ;
    result_x = calloc (Ax_size * sizeof (int64_t)) ;
    memcpy (result_p, Ap, n * sizeof (int64_t)) ;
    memcpy (result_j, Aj, Aj_size * sizeof (int64_t)) ;

    // 2. for ∀s ∈ V
    for (int64_t s = 0; s < n; s++) {
        // 4. S ← empty stack
        size_t sp = 0;

        // Initialize predecessors list P[w] to empty
        // TODO
        // 5. P [w] ← empty queue, ∀w ∈ V
        Pj = malloc (n * sizeof (int64_t)) ;
        Ptail = malloc (n * sizeof (int64_t)) ;
        Phead = Ap ;
        memcpy (Ptail, Ap, n * sizeof (int64_t)) ;

        // Initialize sigma[t], d[t] for all t
        // 6. σ[t] ← 0, ∀t ∈ V , σ[s] ← 1
        // Keeps track of the number of shortest paths for each vertex.
        int64_t *sigma = calloc(n, sizeof(int64_t)) ;
        sigma[s] = 1 ;

        // 7. d[t] ← −1, ∀t ∈ V , d[s] ← 0
        for (size_t t = 0; t < n; t++) {
            d[t] = -1;
        }
        d[s] = 0;

        // Initialize queue and enqueue starting node s
        // 8. Q ← empty queue
        int64_t qh = 0, qt = 0;
        // 9. enqueue(Q, s)
        queue[0] = s;

        // 10. while ¬empty(Q)
        while (qh < qt) {
            // Dequeue v from Q and push onto S
            // 12. v ← dequeue(Q)
            int64_t v = queue[qh++] ;

            // 13. push(S, v)
            S[sp++] = v;

            // TODO
            // traverse all entries in A(v,:)
            for (int64_t p = Ap [v] ; p < Ap [v+1] ; p++)
            {
                int64_t w = Aj [p] ;
                
                // 16. if d[w] < 0
                if (d[w] < 0) {
                    // Update depth and enqueue
                    // 18. enqueue(Q, w)
                    queue[qt++] = w;
                    // 19. d[w] ← d[v] + 1
                    d[w] = d[v] + 1
                }

                // 20. if d[w] = d[v] + 1
                if (d[w] == d[v] + 1) {
                    // Update shortest path count and add predecessor
                    // 22. σ[w] ← σ[w] + σ[v]
                    sigma[w] = sigma[w] + sigma[v]
                    // 23. append(P [w], v)
                    Pj [Ptail [w]++] = v ;
                }

            }       
        }

        // Set dependency score δ[v] ← 0
        // 24. δ[v] ← 0, ∀v ∈ V
        for (size_t v = 0; v < n; v++) {
            d[v] = 0;
        }

        // Process stack S
        // 25. while ¬empty(S)
        while (sp > 0) {
            // 27. w ← pop(S)
            int64_t w = S[--sp];

            // 28. for v ∈ P [w]
            for (int64_t p = Phead [w] ; p < Ptail [w+1] ; p++)
            {
                int64_t v = Pj [p] ;
                
                // Update dependency and centrality values
                // 30. δ[v] ← δ[v] + σ[v] × ( δ[w]/σ[w] + 1)
                double centrality = sigma[v] * ((delta[w] / sigma[w]) + 1);
                delta[v] += centrality;

                // 31. result [(v, w)] ← result [(v, w)] + σ[v] × ( δ[w]/σ[w] + 1)
                size_t w_i = 0;
                for (size_t i = result_p[v]; i < result_p[v + 1]; i++) {
                    if (result_j[i] == w) {
                        w_i = i - result_p[v];
                        break ;
                    }
                }

                result_x[result_p[v] + w_i] += centrality;

            }
        }
    }

    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("LG_check_bfs bfs   time: %g sec\n", tt) ;
        tt = LAGraph_WallClockTime ( ) ;
    }

    //--------------------------------------------------------------------------
    // repack the A matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_Matrix_pack_CSR (G->A,
        &Ap, &Aj, &Ax, Ap_size, Aj_size, Ax_size, iso, jumbled, NULL)) ;
    #endif

    (*C) = result ;

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
}
