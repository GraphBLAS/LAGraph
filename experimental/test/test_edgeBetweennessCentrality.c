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
    LAGraph_Free ((void **) &queue, NULL) ;         \
    LAGraph_Free ((void **) &level_check, NULL) ;   \
    LAGraph_Free ((void **) &level_in, NULL) ;      \
    LAGraph_Free ((void **) &parent_in, NULL) ;     \
    LAGraph_Free ((void **) &visited, NULL) ;       \
    LAGraph_Free ((void **) &neighbors, NULL) ;     \
    GrB_free (&Row) ;                               \
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

    // A centrality matrix initialized to 0 for all edges.
    GrB_Matrix CB = NULL ;

    // Keeps track of the number of shortest paths.
    GrB_Vector sigma = NULL ;

    // Holds the distances (depth levels) from the source vertex.
    GrB_Vector d = NULL ;

    // Stores dependency scores for each vertex.
    GrB_Vector delta = NULL ;

    GrB_Index *S = NULL, *P = NULL ;

    LG_TRY (LAGraph_Malloc ((void **) &queue, n, sizeof (int64_t), msg)) ;

    //--------------------------------------------------------------------------
    // unpack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    bool iso, jumbled ;
    GRB_TRY (GxB_Matrix_unpack_CSR (G->A,
        &Ap, &Aj, &Ax, &Ap_size, &Aj_size, &Ax_size, &iso, &jumbled, NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // compute the level of each node
    //--------------------------------------------------------------------------

    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("LG_check_bfs init  time: %g sec\n", tt) ;
        tt = LAGraph_WallClockTime ( ) ;
    }

    // Initialize centrality matrix CB to 0
    // 1. CB [(v, w)] ← 0, ∀(v, w) ∈ E
    GRB_TRY (GrB_Matrix_new(CB, GrB_FP64, n, n)) ;

    GRB_TRY (GrB_Vector_new(&sigma, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new(&d, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_new(&delta, GrB_FP64, n)) ;

    LG_TRY(LAGraph_Malloc((void **) &S, n, sizeof(GrB_Index), msg)) ;
    LG_TRY(LAGraph_Malloc((void **) &P, n * n, sizeof(GrB_Index), msg)) ;

    GrB_Index queue[n];

    // 2. for ∀s ∈ V
    for (GrB_Index s = 0; s < n; s++) {
        // 4. S ← empty stack
        size_t sp = 0;

        // Initialize predecessors list P[w] to empty
        // 5. P [w] ← empty list, ∀w ∈ V
        memset(P, 0, n * n * sizeof(GrB_Index));

        // Initialize sigma[t], d[t] for all t
        // 6. σ[t] ← 0, ∀t ∈ V , σ[s] ← 1
        GRB_TRY (GrB_Vector_setElement(sigma, 1, s));
        // 7. d[t] ← −1, ∀t ∈ V , d[s] ← 0
        GRB_TRY (GrB_Vector_setElement(d, 0, s));

        // Initialize queue and enqueue starting node s
        // 8. Q ← empty queue
        size_t qp = 0, qlen = 1;
        // 9. enqueue(Q, s)
        queue[0] = s;

        // 10. while ¬empty(Q)
        while (qlen > 0) {
            // Dequeue v from Q and push onto S
            // 12. v ← dequeue(Q)
            GrB_Index v = queue[qp++];
            qlen--;

            // 13. push(S, v)
            S[sp++] = v;

            // Iterate over neighbors of v
            GrB_Vector v_neighbors = NULL;
            GRB_TRY (GrB_Vector_new(&v_neighbors, GrB_BOOL, n));
            GrB_Col_extract(v_neighbors, NULL, NULL, A, GrB_ALL, n, v, GrB_DESC_T0);

            GrB_Index w;
            GrB_Index nvals;
            GrB_Vector_nvals(&nvals, v_neighbors);
            GrB_Index *neighbors = malloc(nvals * sizeof(GrB_Index));
            GrB_Vector_extractTuples_BOOL(neighbors, NULL, &nvals, v_neighbors);

            // 14. for ∀w ∈ neighbors(v)
            for (GrB_Index i = 0; i < nvals; i++) {
                w = neighbors[i];
                int64_t d_w;
                bool d_w_exists = GrB_Vector_extractElement_INT64(&d_w, d, w) == GrB_SUCCESS;

                // 16. if d[w] < 0
                if (!d_w_exists || d_w < 0) {
                    // Update depth and enqueue
                    // 18. enqueue(Q, w)
                    queue[qp + qlen++] = w;
                    // 19. d[w] ← d[v] + 1
                    GrB_Vector_setElement(d, d[v] + 1, w);
                }

                // 20. if d[w] = d[v] + 1
                if (d_w == d[v] + 1) {
                    // Update shortest path count and add predecessor
                    // 22. σ[w] ← σ[w] + σ[v]
                    GrB_Vector_setElement(sigma, sigma[v] + sigma[w], w);
                    // 23. append(P [w], v)
                    P[w * n + v] = 1;
                }
            }


            GrB_Vector_free(&v_neighbors);
            free(neighbors);
        }

        // Set dependency score δ[v] ← 0
        // 24. δ[v] ← 0, ∀v ∈ V
        GrB_Vector_clear(delta);

        // Process stack S
        // 25. while ¬empty(S)
        while (sp > 0) {
            // 27. w ← pop(S)
            GrB_Index w = S[--sp];

            // 28. for v ∈ P [w]
            for (GrB_Index v = 0; v < n; v++) {
                if (P[w * n + v]) {
                    // Update dependency and centrality values
                    double sigma_v, sigma_w, delta_w;
                    GrB_Vector_extractElement(&sigma_v, sigma, v);
                    GrB_Vector_extractElement(&sigma_w, sigma, w);
                    GrB_Vector_extractElement(&delta_w, delta, w);

                    double contribution = sigma_v * ((delta_w / sigma_w) + 1);
                    // 30. δ[v] ← δ[v] + σ[v] × ( δ[w]/σ[w] + 1)
                    delta[v] += contribution;
                    // 31. CB [(v, w)] ← CB [(v, w)] + σ[v] × ( δ[w]/σ[w] + 1)
                    CB[v * n + w] += contribution;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // TODO: get rid of this?
    // for (x = P_head(i); x != 1; x = P_next(x))
    // {
    //     for (p = Ap[i]; p < Ap[i + 1]; p++)
    // }
    //--------------------------------------------------------------------------


    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("LG_check_bfs bfs   time: %g sec\n", tt) ;
        tt = LAGraph_WallClockTime ( ) ;
    }

    //--------------------------------------------------------------------------
    // repack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_Matrix_pack_CSR (G->A,
        &Ap, &Aj, &Ax, Ap_size, Aj_size, Ax_size, iso, jumbled, NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;

    if (print_timings)
    {
        tt = LAGraph_WallClockTime ( ) - tt ;
        printf ("test_edgeBetweennessCentrality check time: %g sec\n", tt) ;
    }
    return (GrB_SUCCESS) ;
}
