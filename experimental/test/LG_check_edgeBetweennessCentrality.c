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
    LAGraph_Free ((void **) &d, NULL) ;             \
    LAGraph_Free ((void **) &delta, NULL) ;         \
    LAGraph_Free ((void **) &S, NULL) ;             \
    LAGraph_Free ((void **) &sigma, NULL) ;         \
    LAGraph_Free ((void **) &Pj, NULL) ;            \
    LAGraph_Free ((void **) &Ptail, NULL) ;         \
}

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
    LAGraph_Free ((void **) &Ap, NULL) ;            \
    LAGraph_Free ((void **) &Aj, NULL) ;            \
    LAGraph_Free ((void **) &Ax, NULL) ;            \
}

#include "LG_internal.h"
#include <LAGraphX.h>

//------------------------------------------------------------------------------
// test the results from a Edge Betweenness Centrality
//------------------------------------------------------------------------------

int LG_check_edgeBetweennessCentrality
(
    // output
    GrB_Matrix *C,      // centrality matrix
    // input
    LAGraph_Graph G,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------

    double tt = LAGraph_WallClockTime ( ) ;

    GrB_Info info;

    double* result ; 

    // Holds the distances (depth levels) from the source vertex.
    int64_t *d = NULL ;

    // Stores dependency scores for each vertex.
    double *delta = NULL ;

    // Stack used for backtracking phase
    int64_t *S = NULL ;

    // Queue used for BFS phase
    int64_t *queue = NULL ;

    GrB_Index *Pj = NULL ;
    GrB_Index *Ptail = NULL ;
    GrB_Index *Phead = NULL ;

    double *sigma = NULL ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Index *Ap = NULL, *Aj = NULL, *neighbors = NULL ;
    void *Ax = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, nvals ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
    bool print_timings = (n >= 2000) ;

    GrB_Matrix A = G->A ;

    GrB_Matrix AT ;
    // if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
    //     G->is_symmetric_structure == LAGraph_TRUE)
    // {
    //     // A and A' have the same structure
    //     AT = A ;
    // }
    // else
    // {
    //     // A and A' differ
    //     AT = G->AT ;
    //     LG_ASSERT_MSG (AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required") ;
    // }

    // better: as basic algo
    /*
        GrB_Matrix A_temp = NULL
        if nself_edges is unknown
            compute it
        if any self edges
            copy G->A into A_temp
            remove self edges from A_temp
            A = A_temp
        now A has no self edges
        when done: free A_temp
    */

    // hack:
    // G->nself_edges = LAGRAPH_UNKNOWN ; <=== overkill
    LG_TRY (LAGraph_DeleteSelfEdges (G, msg)) ;

    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    
    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    LG_TRY(LAGraph_Malloc((void **)&d, n, sizeof(int64_t), msg));

    LG_TRY(LAGraph_Calloc((void **)&delta, n, sizeof(double), msg));

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
    // TODO make this a copy of A except with 1 = 0
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
    bool iso ; 
    GRB_TRY (GxB_Matrix_unpack_CSR (A,
        &Ap, &Aj, &Ax, &Ap_size, &Aj_size, &Ax_size, &iso, NULL, NULL)) ;
    #endif

    Phead = Ap ;

    //--------------------------------------------------------------------------

    LG_TRY(LAGraph_Malloc((void **)&Pj, nvals, sizeof(GrB_Index), msg));
    LG_TRY(LAGraph_Malloc((void **)&Ptail, n, sizeof(GrB_Index), msg)); // might need to be + 1

    LAGraph_Calloc ((void **) &sigma, n, sizeof (double), msg) ;

    // 2. for ∀s ∈ V
    for (int64_t s = 0; s < n; s++) {
        // 4. S ← empty stack
        size_t sp = 0;

        // Initialize predecessors list P[w] to empty
        // TODO
        // 5. P [w] ← empty queue, ∀w ∈ V
        memcpy (Ptail, Ap, n * sizeof (GrB_Index)) ;

        // Initialize sigma[t], d[t] for all t
        // 6. σ[t] ← 0, ∀t ∈ V , σ[s] ← 1
        // Keeps track of the number of shortest paths for each vertex.
        for (int64_t i = 0; i < n; i++) {
            sigma [i] = 0 ;
        }
        sigma [s] = 1 ;

        // 7. d[t] ← −1, ∀t ∈ V , d[s] ← 0
        for (size_t t = 0; t < n; t++) {
            d [t] = -1;
        }
        d [s] = 0;

        // Initialize queue and enqueue starting node s
        // 8. Q ← empty queue
        int64_t qh = 0, qt = 0;
        // 9. enqueue(Q, s)
        queue [qt++] = s;

        // 10. while ¬empty(Q)
        while (qh < qt) {
            // Dequeue v from Q and push onto S
            // 12. v ← dequeue(Q)
            int64_t v = queue [qh++] ;

            // 13. push(S, v)
            S [sp++] = v;

            // TODO
            // traverse all entries in A(v,:)
            for (int64_t p = Ap [v] ; p < Ap [v+1] ; p++)
            {
                int64_t w = Aj [p] ;
                
                // 16. if d[w] < 0
                if (d [w] < 0) {
                    // Update depth and enqueue
                    // 18. enqueue(Q, w)
                    queue [qt++] = w ;
                    // 19. d[w] ← d[v] + 1
                    d [w] = d [v] + 1 ;
                }

                // 20. if d[w] = d[v] + 1
                if (d [w] == d [v] + 1) {
                    // Update shortest path count and add predecessor
                    // 22. σ[w] ← σ[w] + σ[v]
                    sigma [w] = sigma [w] + sigma [v] ;
                    // 23. append(P [w], v)
                    Pj [Ptail [w]++] = v ;
                }

            }       
        }

        // Set dependency score δ[v] ← 0
        // 24. δ[v] ← 0, ∀v ∈ V
        for (size_t v = 0; v < n; v++) {
            delta [v] = 0 ;
        }

        // Process stack S
        // 25. while ¬empty(S)
        while (sp > 0) {
            // 27. w ← pop(S)
            int64_t w = S [--sp] ;

            // 28. for v ∈ P [w]
            for (int64_t p = Phead [w] ; p < Ptail [w] ; p++)
            {
                int64_t v = Pj [p] ;
                
                // Update dependency and centrality values
                // 30. δ[v] ← δ[v] + σ[v] × ( δ[w]/σ[w] + 1)
                printf("%g = %g * (%g/%g + 1)\n", sigma [v] * ((delta [w] / sigma [w]) + 1), sigma [v], delta [w], sigma [w]) ;

                // if (v == w) { printf ("Ack!!\n") ; fflush (stdout) ; abort ( ) ; }
                if (v == w) { 
                    printf ("Ack!!\n") ; 
                    goto flag;
                }

                double centrality = sigma [v] * ((delta [w] / sigma [w]) + 1) ;
                delta [v] += centrality ;

                // 31. result [(v, w)] ← result [(v, w)] + σ[v] × ( δ[w]/σ[w] + 1)
                result [INDEX (v,w)] += centrality;

            }
        }

    }

    flag:
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
        &Ap, &Aj, &Ax, Ap_size, Aj_size, Ax_size, iso, NULL, NULL)) ;
    #endif

    printf("result: \n") ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        printf ("row: %ld\n", i) ;
        for (int64_t j = 0 ; j < n ; j++)
        {
            int64_t p = INDEX (i,j) ;
            double aij = result [p] ;
            printf("  C(%ld,%ld) = %g\n", i, j, aij) ;
            // numerical value of A(i,j)
        } 
        printf("\n") ; 
    }

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

    GxB_print(C_temp, GxB_COMPLETE) ;

    // GrB_TRY (GrB_select(*C_temp, A, NULL, NULL, A, NULL, NULL)) ;

    *C = C_temp;

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
