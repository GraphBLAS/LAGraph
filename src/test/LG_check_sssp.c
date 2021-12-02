//------------------------------------------------------------------------------
// LAGraph/src/test/LG_check_sssp: stand-alone test for SSSP
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LG_test.h"

// TODO: handle all data types, not just int32

typedef int64_t LG_key_t ;
typedef struct
{
    int64_t name ;
    LG_key_t key ;
}
LG_Element ;
#include "LG_heap.h"

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK                           \
{                                                   \
    LAGraph_Free ((void **) &Heap) ;                \
    LAGraph_Free ((void **) &Iheap) ;               \
    LAGraph_Free ((void **) &distance) ;            \
    LAGraph_Free ((void **) &parent) ;              \
    LAGraph_Free ((void **) &path_length_in) ;      \
    LAGraph_Free ((void **) &neighbor_weights) ;    \
    LAGraph_Free ((void **) &neighbors) ;           \
    GrB_free (&Row) ;                               \
}

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL                    \
{                                           \
    LAGraph_FREE_WORK ;                     \
    LAGraph_Free ((void **) &Ap) ;          \
    LAGraph_Free ((void **) &Aj) ;          \
    LAGraph_Free ((void **) &Ax) ;          \
}

//------------------------------------------------------------------------------
// test the results from SSSP
//------------------------------------------------------------------------------

// Because this method does on GxB_unpack on G->A, it should not be used in a
// brutal memory test, unless the caller is prepared to reconstruct G->A
// when the brutal test causes this method to return early.

int LG_check_sssp
(
    // input
    GrB_Vector Path_Length,     // Path_Length(i) is the length of the
                                // shortest path from src to node i.
    LAGraph_Graph G,            // all edge weights must be > 0
    GrB_Index src,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    double tic [2], tt ;
    GrB_Vector Row = NULL ;
    GrB_Index *Ap = NULL, *Aj = NULL, *neighbors = NULL ;
    int32_t *Ax = NULL, *neighbor_weights = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, ncols ;
    int64_t *queue = NULL, *path_length_in = NULL, *Iheap = NULL,
        *distance = NULL, *parent = NULL ;
    LG_Element *Heap = NULL ;

    LAGraph_Tic (tic, msg) ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_CHECK (n != ncols, -1001, "G->A must be square") ;
    LG_CHECK (G->A_type != GrB_INT32, -1001, "G->A must be int32") ;
    bool print_timings = (n >= 2000) ;

    //--------------------------------------------------------------------------
    // get the contents of the Path_Length vector
    //--------------------------------------------------------------------------

    path_length_in = LAGraph_Malloc (n, sizeof (int64_t)) ;
    LG_CHECK (path_length_in == NULL, -1003, "out of memory") ;
    LG_CHECK (LG_check_vector (path_length_in, Path_Length, n, INT32_MAX),
        -1004, "invalid Path_Length") ;

    //--------------------------------------------------------------------------
    // unpack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    bool iso = false ;
    #if LG_SUITESPARSE
    bool jumbled ;
    GrB_TRY (GxB_Matrix_unpack_CSR (G->A,
        &Ap, &Aj, (void **) &Ax, &Ap_size, &Aj_size, &Ax_size, &iso, &jumbled,
        NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // compute the SSSP of the graph, via Dijskstra's algorithm
    //--------------------------------------------------------------------------

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_sssp init  time: %g sec\n", tt) ;
        LAGraph_Tic (tic, msg) ;
    }

    // initializations
    distance = LAGraph_Malloc (n, sizeof (int64_t)) ;
    // parent = LAGraph_Malloc (n, sizeof (int64_t)) ;
    LG_CHECK (distance == NULL, -1003, "out of memory") ;
    // LG_CHECK (parent == NULL, -1003, "out of memory") ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        distance [i] = INT32_MAX ;
        // parent [i] = -1 ;
    }
    distance [src] = 0 ;
    // parent [src] = src ;

    #if !LG_SUITESPARSE
    GrB_TRY (GrB_Vector_new (&Row, GrB_INT32, n)) ;
    neighbors = LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    neighbor_weights = LAGraph_Malloc (n, sizeof (int32_t)) ;
    LG_CHECK (neighbors == NULL, -1003, "out of memory") ;
    LG_CHECK (neighbor_weights == NULL, -1003, "out of memory") ;
    #endif

    // place all nodes in the heap (already in heap order)
    Heap = LAGraph_Malloc ((n+1), sizeof (LG_Element)) ;
    Iheap = LAGraph_Malloc (n, sizeof (int64_t)) ;
    LG_CHECK (Heap == NULL || Iheap == NULL, -1003, "out of memory") ;
    Heap [1].key = 0 ;
    Heap [1].name = src ;
    Iheap [src] = 1 ;
    int64_t p = 2 ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        if (i != src)
        {
            Heap [p].key = INT32_MAX ;
            Heap [p].name = i ;
            Iheap [i] = p ;
            p++ ;
        }
    }
    int64_t nheap = n ;
    LG_CHECK (LG_heap_check (Heap, Iheap, n, nheap), -1011, "invalid heap") ;

    while (nheap > 0)
    {
        // extract the min element u from the top of the heap
        LG_Element e = Heap [1] ;
        int64_t u = e.name ;
        // printf ("\nGot node %ld\n", u) ;

        int64_t u_distance = e.key ;
        ASSERT (distance [u] == u_distance) ;
        LG_heap_delete (1, Heap, Iheap, n, &nheap) ;
        ASSERT (Iheap [u] == 0) ;

        // printf ("\nafter delete\n") ;
        if (n < 200)
        {
            LG_CHECK (LG_heap_check (Heap, Iheap, n, nheap), -1013,
                "invalid heap") ;
        }

        if (u_distance == INT32_MAX)
        {
            // node u is not reachable, so no other nodes in the queue
            // are reachable either.  All work is done.
            break ;
        }

        #if LG_SUITESPARSE
        // directly access the indices of entries in A(u,:)
        GrB_Index degree = Ap [u+1] - Ap [u] ;
        GrB_Index *node_u_adjacency_list = Aj + Ap [u] ;
        int32_t *weights = Ax + (iso ? 0 : Ap [u]) ;
        #else
        // extract the indices of entries in A(u,:)
        GrB_Index degree = n ;
        GrB_TRY (GrB_Col_extract (Row, NULL, NULL, G->A, GrB_ALL, n, u,
            GrB_DESC_T0)) ;
        GrB_TRY (GrB_Vector_extractTuples_INT32 (neighbors, neighbor_weights,
            &degree, Row)) ;
        GrB_Index *node_u_adjacency_list = neighbors ;
        int32_t *weights = neighbor_weights ;
        #endif

        #if 0
        printf ("Ap %ld %ld\n", Ap [u], Ap [u+1]) ;
        // check all entries in A(u,:)
        for (int64_t k = 0 ; k < degree ; k++)
        {
            // consider edge (u,v) and its weight w
            int64_t v = node_u_adjacency_list [k] ;
            int64_t w = (int64_t) (weights [iso ? 0 : k]) ;
            printf ("consider edge u %ld (in heap: %ld)\n", u, Iheap [u]) ;
            printf ("consider node v %ld (in heap: %ld)\n", v, Iheap [v]) ;
            printf ("weight w %ld\n", w) ;
        }
        printf ("\n") ;
        #endif

        // traverse all entries in A(u,:)
        for (int64_t k = 0 ; k < degree ; k++)
        {
            // consider edge (u,v) and its weight w
            int64_t v = node_u_adjacency_list [k] ;
            if (Iheap [v] == 0) continue ;  // node v already in SSSP tree
            int64_t w = (int64_t) (weights [iso ? 0 : k]) ;
            // printf ("consider edge (%ld,%ld) weight %ld\n", u, v, w) ;

            LG_CHECK (w <= 0, -1008, "invalid graph (weights must be > 0)") ;
            int64_t new_distance = u_distance + w ;
            if (distance [v] > new_distance)
            {
                // reduce the key of node v
                // printf ("decreased key of node %ld from %ld to %ld\n",
                //     v, distance [v], new_distance) ;
                distance [v] = new_distance ;
                // parent [v] = u ;
                int64_t p = Iheap [v] ;
                LG_CHECK (Heap [p].name != v, -2000, "huh?") ;
                LG_heap_decrease_key (p, new_distance, Heap, Iheap, n, nheap) ;
            }
        }

        if (n < 200)
        {
            LG_CHECK (LG_heap_check (Heap, Iheap, n, nheap), -1014,
                "invalid heap") ;
        }

    }

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_sssp time: %g sec\n", tt) ;
        LAGraph_Tic (tic, msg) ;
    }

    //--------------------------------------------------------------------------
    // repack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    GrB_TRY (GxB_Matrix_pack_CSR (G->A,
        &Ap, &Aj, (void **) &Ax, Ap_size, Aj_size, Ax_size, iso, jumbled,
        NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // check the distance of each node
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; i < n ; i++)
    {
        bool ok = (path_length_in [i] == distance [i]) ;
        LG_CHECK (!ok, -1004, "invalid path length") ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_sssp check time: %g sec\n", tt) ;
    }
    return (0) ;
}

