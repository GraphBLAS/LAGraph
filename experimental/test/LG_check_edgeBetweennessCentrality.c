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

    // A temporary result centrality matrix initialized to 0 for all edges.
    GrB_Matrix result = NULL ;

    // Keeps track of the number of shortest paths.
    GrB_Vector sigma = NULL ;

    // Holds the distances (depth levels) from the source vertex.
    GrB_Vector d = NULL ;

    // Stores dependency scores for each vertex.
    GrB_Vector delta = NULL ;

    GrB_Index *S = NULL, *P = NULL ;

    LG_TRY (LAGraph_Malloc ((void **) &queue, n, sizeof (int64_t), msg)) ;

    //--------------------------------------------------------------------------
    // unpack the A matrix in CSR form for SuiteSparse:GraphBLAS
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

    // Initialize centrality matrix result to 0
    // 1. result [(v, w)] ← 0, ∀(v, w) ∈ E
    GRB_TRY (GrB_Matrix_new(result, GrB_FP64, n, n, A)) ;
    GRB_TRY (GrB_assign(result, A, null, 0, GrB_ALL, n, Grb_ALL, n, GrB_DESC_S)) ;

    //--------------------------------------------------------------------------
    // unpack the centrality matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    bool Ciso, Cjumbled ;
    GRB_TRY (GxB_Matrix_unpack_CSR (result,
        &Cp, &Cj, &Cx, &Cp_size, &Cj_size, &Cx_size, &Ciso, &Cjumbled, NULL)) ;
    #endif


    GRB_TRY (GrB_Vector_new(&sigma, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new(&d, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_new(&delta, GrB_FP64, n)) ;

    LG_TRY(LAGraph_Malloc((void **) &S, n, sizeof(GrB_Index), msg)) ;
    LG_TRY(LAGraph_Malloc((void **) &P, n * n, sizeof(GrB_Index), msg)) ;

    #if !LAGRAPH_SUITESPARSE
    GRB_TRY (GrB_Vector_new (&Row, GrB_BOOL, n)) ;
    LG_TRY (LAGraph_Malloc ((void **) &neighbors, n, sizeof (GrB_Index), msg)) ;
    #endif

    // 2. for ∀s ∈ V
    for (int64_t s = 0; s < n; s++) {
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
        int64_t qp = 0, qlen = 1;
        // 9. enqueue(Q, s)
        queue[0] = s;

        // 10. while ¬empty(Q)
        while (qlen > 0) {
            // Dequeue v from Q and push onto S
            // 12. v ← dequeue(Q)
            int64_t v = queue [qp++] ;
            qlen--;

            // 13. push(S, v)
            S[sp++] = v;

            #if LAGRAPH_SUITESPARSE
            // directly access the indices of entries in A(v,:)
            GrB_Index degree = Ap [v+1] - Ap [v] ;
            GrB_Index *node_u_adjacency_list = Aj + Ap [v] ;
            #else
            // extract the indices of entries in A(v,:)
            GrB_Index degree = n ;
            GRB_TRY (GrB_Col_extract (Row, NULL, NULL, G->A, GrB_ALL, n, v,
                GrB_DESC_T0)) ;
            GRB_TRY (GrB_Vector_extractTuples_BOOL (neighbors, NULL, &degree, Row));
            GrB_Index *node_v_adjacency_list = neighbors ;
            #endif

            // traverse all entries in A(v,:)
            for (int64_t k = 0 ; k < degree ; k++)
            {
                // consider edge (v,w)
                int64_t w = node_v_adjacency_list [k] ;
                int64_t d_w;
                bool d_w_exists = GrB_Vector_extractElement_INT64(&d_w, d, w) == GrB_SUCCESS;

                // 16. if d[w] < 0
                if (!d_w_exists || d_w < 0) {
                    // Update depth and enqueue
                    // 18. enqueue(Q, w)
                    queue[qp + qlen++] = w;
                    // 19. d[w] ← d[v] + 1
                    GRB_TRY (GrB_Vector_setElement(d, d[v] + 1, w)) ;
                }

                // 20. if d[w] = d[v] + 1
                if (d_w == d[v] + 1) {
                    // Update shortest path count and add predecessor
                    // 22. σ[w] ← σ[w] + σ[v]
                    GRB_TRY (GrB_Vector_setElement(sigma, sigma[v] + sigma[w], w)) ;
                    // 23. append(P [w], v)
                    P[w * n + v] = 1;
                }
            }
        }

        // Set dependency score δ[v] ← 0
        // 24. δ[v] ← 0, ∀v ∈ V
        GRB_TRY (GrB_Vector_clear(delta));

        // Process stack S
        // 25. while ¬empty(S)
        while (sp > 0) {
            // 27. w ← pop(S)
            int64_t w = S[--sp];

            // 28. for v ∈ P [w]
            for (int64_t v = 0; v < n; v++) {
                if (P[w * n + v]) {
                    // Update dependency and centrality values
                    double sigma_v, sigma_w, delta_w;
                    GRB_TRY (GrB_Vector_extractElement(&sigma_v, sigma, v));
                    GRB_TRY (GrB_Vector_extractElement(&sigma_w, sigma, w));
                    GRB_TRY (GrB_Vector_extractElement(&delta_w, delta, w));

                    double centrality = sigma_v * ((delta_w / sigma_w) + 1);
                    // 30. δ[v] ← δ[v] + σ[v] × ( δ[w]/σ[w] + 1)
                    delta[v] += centrality;
                    // 31. result [(v, w)] ← result [(v, w)] + σ[v] × ( δ[w]/σ[w] + 1)
                    result[v * n + w] += centrality;

                    // TODO: maybe get rid of this
                    // int x;
                    // GrB_extract(&x, c, v, w) ; 
                }
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

     //--------------------------------------------------------------------------
    // repack the centrality matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_Matrix_pack_CSR (result,
        &Cp, &Cj, &Cx, Cp_size, Cj_size, Cx_size, Ciso, Cjumbled, NULL)) ;
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
