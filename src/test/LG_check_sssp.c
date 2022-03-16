//------------------------------------------------------------------------------
// LAGraph/src/test/LG_check_sssp: stand-alone test for SSSP
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include "LG_test.h"

// All computations are done in double precision

typedef double LG_key_t ;
typedef struct
{
    int64_t name ;
    LG_key_t key ;
}
LG_Element ;
#include "LG_heap.h"

#undef  LG_FREE_WORK
#define LG_FREE_WORK                                    \
{                                                       \
    LAGraph_Free ((void **) &Heap, NULL) ;              \
    LAGraph_Free ((void **) &Iheap, NULL) ;             \
    LAGraph_Free ((void **) &distance, NULL) ;          \
    LAGraph_Free ((void **) &parent, NULL) ;            \
    LAGraph_Free ((void **) &path_length_in, NULL) ;    \
    LAGraph_Free ((void **) &neighbor_weights, NULL) ;  \
    LAGraph_Free ((void **) &neighbors, NULL) ;         \
    GrB_free (&Row) ;                                   \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                                     \
{                                                       \
    LG_FREE_WORK ;                                      \
    LAGraph_Free ((void **) &Ap, NULL) ;                \
    LAGraph_Free ((void **) &Aj, NULL) ;                \
    LAGraph_Free ((void **) &Ax, NULL) ;                \
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
    GrB_Index Ap_size, Aj_size, Ax_size, n, ncols ;
    int64_t *queue = NULL, *Iheap = NULL, *parent = NULL ;
    LG_Element *Heap = NULL ;

    double *path_length_in = NULL, *distance = NULL, *neighbor_weights = NULL ;
    void *Ax = NULL ;

    LAGraph_Tic (tic, msg) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    GrB_Type etype ;
    char atype_name [LAGRAPH_MAX_NAME_LEN] ;
    LG_TRY (LAGraph_Matrix_TypeName (atype_name, G->A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&etype, atype_name, msg)) ;
    int etypecode = 0 ;
    size_t etypesize = 0 ;
    double etypeinf = INFINITY ;

    if (etype == GrB_INT32)
    {
        etypecode = 0 ;
        etypesize = sizeof (int32_t) ;
        etypeinf  = INT32_MAX ;
    }
    else if (etype == GrB_INT64)
    {
        etypecode = 1 ;
        etypesize = sizeof (int64_t) ;
        etypeinf  = INT64_MAX ;
    }
    else if (etype == GrB_UINT32)
    {
        etypecode = 2 ;
        etypesize = sizeof (uint32_t) ;
        etypeinf  = UINT32_MAX ;
    }
    else if (etype == GrB_UINT64)
    {
        etypecode = 3 ;
        etypesize = sizeof (uint64_t) ;
        etypeinf  = UINT64_MAX ;
    }
    else if (etype == GrB_FP32)
    {
        etypecode = 4 ;
        etypesize = sizeof (float) ;
        etypeinf  = INFINITY ;
    }
    else if (etype == GrB_FP64)
    {
        etypecode = 5 ;
        etypesize = sizeof (double) ;
        etypeinf  = INFINITY ;
    }

    LG_ASSERT_MSG (etypesize > 0, GrB_NOT_IMPLEMENTED, "type not supported") ;

    bool print_timings = (n >= 2000) ;

    //--------------------------------------------------------------------------
    // get the contents of the Path_Length vector
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Malloc ((void **) &path_length_in, n, sizeof (double), msg)) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        double t ;
        path_length_in [i] = INFINITY ;
        int info = GrB_Vector_extractElement_FP64 (&t, Path_Length, i) ;
        if (info == GrB_SUCCESS)
        {
            path_length_in [i] = t ;
        }
    }

    //--------------------------------------------------------------------------
    // unpack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    bool iso = false ;
    #if LAGRAPH_SUITESPARSE
    bool jumbled ;
    GRB_TRY (GxB_Matrix_unpack_CSR (G->A,
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
    LG_TRY (LAGraph_Malloc ((void **) &distance, n, sizeof (double), msg)) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        distance [i] = INFINITY ;
    }
    distance [src] = 0 ;

    #if !LAGRAPH_SUITESPARSE
    GRB_TRY (GrB_Vector_new (&Row, GrB_FP64, n)) ;
    LG_TRY (LAGraph_Malloc ((void **) &neighbors, n, sizeof (GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &neighbor_weights, n, sizeof (double), msg)) ;
    LG_ASSERT (neighbors != NULL, GrB_OUT_OF_MEMORY) ;
    LG_ASSERT (neighbor_weights != NULL, GrB_OUT_OF_MEMORY) ;
    #endif

    // place all nodes in the heap (already in heap order)
    LG_TRY (LAGraph_Malloc ((void **) &Heap, (n+1), sizeof (LG_Element), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &Iheap, n, sizeof (int64_t), msg)) ;
    LG_ASSERT (Heap != NULL && Iheap != NULL, GrB_OUT_OF_MEMORY) ;
    Heap [1].key = 0 ;
    Heap [1].name = src ;
    Iheap [src] = 1 ;
    int64_t p = 2 ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        if (i != src)
        {
            Heap [p].key = INFINITY ;
            Heap [p].name = i ;
            Iheap [i] = p ;
            p++ ;
        }
    }
    int64_t nheap = n ;
    LG_ASSERT_MSG (LG_heap_check (Heap, Iheap, n, nheap) == 0, -2000,
        "invalid heap") ;

    while (nheap > 0)
    {
        // extract the min element u from the top of the heap
        LG_Element e = Heap [1] ;
        int64_t u = e.name ;
        // printf ("\nGot node %ld\n", u) ;

        double u_distance = e.key ;
        ASSERT (distance [u] == u_distance) ;
        LG_heap_delete (1, Heap, Iheap, n, &nheap) ;
        ASSERT (Iheap [u] == 0) ;

        // printf ("\nafter delete\n") ;
        if (n < 200)
        {
            LG_ASSERT_MSG (LG_heap_check (Heap, Iheap, n, nheap) == 0, -2000,
                "invalid heap") ;
        }

        if (u_distance == INFINITY)
        {
            // node u is not reachable, so no other nodes in the queue
            // are reachable either.  All work is done.
            break ;
        }

        #if LAGRAPH_SUITESPARSE
        // directly access the indices of entries in A(u,:)
        GrB_Index degree = Ap [u+1] - Ap [u] ;
        GrB_Index *node_u_adjacency_list = Aj + Ap [u] ;
        void *weights = Ax + ((iso ? 0 : Ap [u]) * etypesize) ;
        #else
        // extract the indices of entries in A(u,:)
        GrB_Index degree = n ;
        GRB_TRY (GrB_Col_extract (Row, NULL, NULL, G->A, GrB_ALL, n, u,
            GrB_DESC_T0)) ;
        GRB_TRY (GrB_Vector_extractTuples_FP64 (neighbors, neighbor_weights,
            &degree, Row)) ;
        GrB_Index *node_u_adjacency_list = neighbors ;
        double *weights = neighbor_weights ;
        #endif

        // traverse all entries in A(u,:)
        for (int64_t k = 0 ; k < degree ; k++)
        {
            // consider edge (u,v) and its weight w
            int64_t v = node_u_adjacency_list [k] ;
            if (Iheap [v] == 0) continue ;  // node v already in SSSP tree
            double w ;
            #if LAGRAPH_SUITESPARSE
            switch (etypecode)
            {
                default:
                case 0: w = (( int32_t *) weights) [iso ? 0 : k] ; break ;
                case 1: w = (( int64_t *) weights) [iso ? 0 : k] ; break ;
                case 2: w = ((uint32_t *) weights) [iso ? 0 : k] ; break ;
                case 3: w = ((uint64_t *) weights) [iso ? 0 : k] ; break ;
                case 4: w = (( float   *) weights) [iso ? 0 : k] ; break ;
                case 5: w = (( double  *) weights) [iso ? 0 : k] ; break ;
            }
            #else
            w = weights [iso ? 0 : k] ;
            #endif
//          printf ("consider edge (%d,%d) weight %g\n", (int) u, (int) v, w) ;

            LG_ASSERT_MSG (w > 0, -2002, "invalid graph (weights must be > 0)");
            double new_distance = u_distance + w ;
            if (distance [v] > new_distance)
            {
                // reduce the key of node v
//              printf ("decreased key of node %d from %g to %g\n",
//                  (int) v, distance [v], new_distance) ;
                distance [v] = new_distance ;
                // parent [v] = u ;
                int64_t p = Iheap [v] ;
                LG_ASSERT_MSG (Heap [p].name == v, -2000, "invalid heap") ;
                LG_heap_decrease_key (p, new_distance, Heap, Iheap, n, nheap) ;
            }
        }

        if (n < 200)
        {
            LG_ASSERT_MSG (LG_heap_check (Heap, Iheap, n, nheap) == 0, -2000,
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

    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GxB_Matrix_pack_CSR (G->A,
        &Ap, &Aj, (void **) &Ax, Ap_size, Aj_size, Ax_size, iso, jumbled,
        NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // check the distance of each node
    //--------------------------------------------------------------------------

    for (int64_t i = 0 ; i < n ; i++)
    {
        bool ok = true ;
        double err = 0 ;
        if (isinf (distance [i]))
        {
            ok = (path_length_in [i] == etypeinf || isinf (path_length_in [i]));
        }
        else
        {
            err = fabs (path_length_in [i] - distance [i]) ;
            if (err > 0) err = err / LAGRAPH_MAX (path_length_in [i], distance [i]) ;
            ok = (err < 1e-5) ;
        }
//      printf ("%d: %g %g err %g ok: %d\n", (int) i,
//          path_length_in [i], distance [i], err, ok) ;
        LG_ASSERT_MSG (ok, -2001, "invalid path length") ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_sssp check time: %g sec\n", tt) ;
    }
    return (GrB_SUCCESS) ;
}

