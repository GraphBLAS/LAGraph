//------------------------------------------------------------------------------
// LAGraph/src/test/LG_check_cc: stand-alone test for CC
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#define LG_FREE_WORK                            \
{                                               \
    LAGraph_Free ((void **) &queue) ;           \
    LAGraph_Free ((void **) &component_in) ;    \
    LAGraph_Free ((void **) &visited) ;         \
    LAGraph_Free ((void **) &neighbors) ;       \
    GrB_free (&Row) ;                           \
}

#define LG_FREE_ALL                             \
{                                               \
    LG_FREE_WORK ;                              \
    LAGraph_Free ((void **) &Ap) ;              \
    LAGraph_Free ((void **) &Aj) ;              \
    LAGraph_Free ((void **) &Ax) ;              \
}

#include "LG_internal.h"
#include "LG_test.h"

// The output of LAGr_ConnectedComponents is a vector Component, where
// Component(i)=s if node i is in the connected compononent whose
// representative node is node s.  If s is a representative, then
// Component(s)=s.  The number of connected components in the graph G is the
// number of representatives.

//------------------------------------------------------------------------------
// test the results from LAGr_ConnectedComponents 
//------------------------------------------------------------------------------

// Because this method does on GxB_unpack on G->A, it should not be used in a
// brutal memory test, unless the caller is prepared to reconstruct G->A
// when the brutal test causes this method to return early.

int LG_check_cc
(
    // input
    GrB_Vector Component,   // Component(i)=s if node is in Component s
    LAGraph_Graph G,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    double tic [2], tt ;
    LAGraph_Tic (tic, msg) ;
    GrB_Vector Row = NULL ;
    GrB_Index *Ap = NULL, *Aj = NULL, *neighbors = NULL ;
    void *Ax = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, ncols ;
    int64_t *queue = NULL, *component_in = NULL ;
    bool *visited = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_ASSERT (Component != NULL, GrB_NULL_POINTER) ;

    LG_ASSERT_MSG ((G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->structure_is_symmetric == LAGraph_TRUE)),
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be known to be symmetric") ;

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    queue = LAGraph_Calloc (n, sizeof (int64_t)) ;
    LG_ASSERT (queue != NULL, GrB_OUT_OF_MEMORY) ;

    //--------------------------------------------------------------------------
    // get the contents of the Component vector
    //--------------------------------------------------------------------------

    component_in = LAGraph_Malloc (n, sizeof (int64_t)) ;
    LG_ASSERT (component_in != NULL, GrB_OUT_OF_MEMORY) ;
    LG_TRY (LG_check_vector (component_in, Component, n, -1)) ;

    //--------------------------------------------------------------------------
    // find the # of connected components, according to Component vector
    //--------------------------------------------------------------------------

    int64_t *count = queue ;        // use queue as workspace
    int64_t ncomp_in = 0 ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        int64_t comp = component_in [i] ; 
        LG_ASSERT (comp >= 0 && comp < n, -2000) ;
        count [comp]++ ;
        if (comp == i)
        {
            // this is the representative of its component
            ncomp_in++ ;
        }
    }
    printf ("# of components: %g\n", (double) ncomp_in) ;

    if (n < 1000)
    {
        for (int64_t i = 0 ; i < n ; i++)
        {
            if (component_in [i] == i)
            {
                printf ("Component %g, size %g\n", (double) i,
                    (double) count [i]) ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // unpack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    bool iso, jumbled ;
    GrB_TRY (GxB_Matrix_unpack_CSR (G->A,
        &Ap, &Aj, &Ax, &Ap_size, &Aj_size, &Ax_size, &iso, &jumbled, NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // find the connected components via repeated BFS
    //--------------------------------------------------------------------------

    LAGraph_Toc (&tt, tic, msg) ;
    printf ("LG_check_cc init  time: %g sec\n", tt) ;
    LAGraph_Tic (tic, msg) ;

    visited = LAGraph_Calloc (n, sizeof (bool)) ;
    LG_ASSERT (visited != NULL, GrB_OUT_OF_MEMORY) ;

    #if !LG_SUITESPARSE
    GrB_TRY (GrB_Vector_new (&Row, GrB_BOOL, n)) ;
    neighbors = LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    LG_ASSERT (neighbors != NULL, GrB_OUT_OF_MEMORY) ;
    #endif

    int64_t ncomp = 0 ;

    for (int64_t src = 0 ; src < n ; src++)
    {
        // skip this node if already visited
        if (visited [src]) continue ;

        // src node is part of a new connected component, comp
        int64_t comp = component_in [src] ;
        ncomp++ ;
        LG_ASSERT_MSG (ncomp <= ncomp_in, -2001, "wrong # of components") ;

        queue [0] = src ;
        int64_t head = 0 ;
        int64_t tail = 1 ;
        visited [src] = true ;      // src is visited

        while (head < tail)
        {
            // dequeue the node at the head of the queue
            int64_t u = queue [head++] ;

            #if LG_SUITESPARSE
            // directly access the indices of entries in A(u,:)
            GrB_Index degree = Ap [u+1] - Ap [u] ;
            GrB_Index *node_u_adjacency_list = Aj + Ap [u] ;
            #else
            // extract the indices of entries in A(u,:)
            GrB_Index degree = n ;
            GrB_TRY (GrB_Col_extract (Row, NULL, NULL, G->A, GrB_ALL, n, u,
                GrB_DESC_T0)) ;
            GrB_TRY (GrB_Vector_extractTuples_BOOL (neighbors, NULL, &degree,
                Row)) ;
            GrB_Index *node_u_adjacency_list = neighbors ;
            #endif

            // traverse all entries in A(u,:)
            for (int64_t k = 0 ; k < degree ; k++)
            {
                // consider edge (u,v)
                int64_t v = node_u_adjacency_list [k] ;
                // ensure v is in the same connected component as the src node
                LG_ASSERT (comp == component_in [u], -2002) ;
                // printf ("    seen: %ld\n", v) ;
                if (!visited [v])
                {
                    // node v is not yet visited; set its level and add to the
                    // end of the queue
                    visited [v] = true ;
                    queue [tail++] = v ;
                }
            }
        }
    }

    LG_ASSERT_MSG (ncomp == ncomp_in, -2001, "wrong # of components") ;

    LAGraph_Toc (&tt, tic, msg) ;
    printf ("LG_check_cc component time: %g sec\n", tt) ;
    LAGraph_Tic (tic, msg) ;

    //--------------------------------------------------------------------------
    // repack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    GrB_TRY (GxB_Matrix_pack_CSR (G->A,
        &Ap, &Aj, &Ax, Ap_size, Aj_size, Ax_size, iso, jumbled, NULL)) ;
    #endif

    LG_FREE_WORK ;

    LAGraph_Toc (&tt, tic, msg) ;
    printf ("LG_check_cc check time: %g sec\n", tt) ;

    return (GrB_SUCCESS) ;
}

