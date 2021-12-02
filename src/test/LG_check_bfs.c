//------------------------------------------------------------------------------
// LAGraph/src/test/LG_check_bfs: stand-alone test for BFS
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#define LAGraph_FREE_WORK                   \
{                                           \
    LAGraph_Free ((void **) &queue) ;       \
    LAGraph_Free ((void **) &level_check) ; \
    LAGraph_Free ((void **) &level_in) ;    \
    LAGraph_Free ((void **) &parent_in) ;   \
    LAGraph_Free ((void **) &visited) ;     \
    LAGraph_Free ((void **) &neighbors) ;   \
    GrB_free (&Row) ;                       \
}

#define LAGraph_FREE_ALL                    \
{                                           \
    LAGraph_FREE_WORK ;                     \
    LAGraph_Free ((void **) &Ap) ;          \
    LAGraph_Free ((void **) &Aj) ;          \
    LAGraph_Free ((void **) &Ax) ;          \
}

#include "LG_internal.h"
#include "LG_test.h"

//------------------------------------------------------------------------------
// test the results from a BFS
//------------------------------------------------------------------------------

// Because this method does on GxB_unpack on G->A, it should not be used in a
// brutal memory test, unless the caller is prepared to reconstruct G->A
// when the brutal test causes this method to return early.

int LG_check_bfs
(
    // input
    GrB_Vector Level,       // optional; may be NULL
    GrB_Vector Parent,      // optional; may be NULL
    LAGraph_Graph G,
    GrB_Index src,
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
    int64_t *queue = NULL, *level_in = NULL, *parent_in = NULL,
        *level_check = NULL ;
    bool *visited = NULL ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_CHECK (n != ncols, -1001, "G->A must be square") ;
    bool print_timings = (n >= 2000) ;

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    queue = LAGraph_Malloc (n, sizeof (int64_t)) ;
    level_check = LAGraph_Malloc (n, sizeof (int64_t)) ; 
    LG_CHECK (queue == NULL || level_check == NULL , -1003, "out of memory") ;

    //--------------------------------------------------------------------------
    // get the contents of the Level and Parent vectors
    //--------------------------------------------------------------------------

    if (Level != NULL)
    {
        level_in = LAGraph_Malloc (n, sizeof (int64_t)) ;
        LG_CHECK (level_in == NULL, -1003, "out of memory") ;
        LG_CHECK (LG_check_vector (level_in, Level, n, -1), -1004,
            "invalid level") ;
    }

    if (Parent != NULL)
    {
        parent_in = LAGraph_Malloc (n, sizeof (int64_t)) ;
        LG_CHECK (parent_in == NULL, -1003, "out of memory") ;
        LG_CHECK (LG_check_vector (parent_in, Parent, n, -1), -1005,
            "invalid parent") ;
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
    // compute the level of each node
    //--------------------------------------------------------------------------

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_bfs init  time: %g sec\n", tt) ;
        LAGraph_Tic (tic, msg) ;
    }

    queue [0] = src ;
    int64_t head = 0 ;
    int64_t tail = 1 ;
    visited = LAGraph_Calloc (n, sizeof (bool)) ;
    LG_CHECK (visited == NULL, -1003, "out of memory") ;
    visited [src] = true ;      // src is visited, and is level 0

    for (int64_t i = 0 ; i < n ; i++)
    {
        level_check [i] = -1 ;
    }
    level_check [src] = 0 ;

    #if !LG_SUITESPARSE
    GrB_TRY (GrB_Vector_new (&Row, GrB_BOOL, n)) ;
    neighbors = LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    LG_CHECK (neighbors == NULL, -1003, "out of memory") ;
    #endif

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
        GrB_TRY (GrB_Vector_extractTuples_BOOL (neighbors, NULL, &degree, Row));
        GrB_Index *node_u_adjacency_list = neighbors ;
        #endif

        // traverse all entries in A(u,:)
        for (int64_t k = 0 ; k < degree ; k++)
        {
            // consider edge (u,v)
            int64_t v = node_u_adjacency_list [k] ;
            if (!visited [v])
            {
                // node v is not yet visited; set its level and add to the
                // end of the queue
                visited [v] = true ;
                level_check [v] = level_check [u] + 1 ;
                queue [tail++] = v ;
            }
        }
    }

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_bfs bfs   time: %g sec\n", tt) ;
        LAGraph_Tic (tic, msg) ;
    }

    //--------------------------------------------------------------------------
    // repack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    GrB_TRY (GxB_Matrix_pack_CSR (G->A,
        &Ap, &Aj, &Ax, Ap_size, Aj_size, Ax_size, iso, jumbled, NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // check the level of each node
    //--------------------------------------------------------------------------

    if (level_in != NULL)
    {
        for (int64_t i = 0 ; i < n ; i++)
        {
            bool ok = (level_in [i] == level_check [i]) ;
            LG_CHECK (!ok, -1004, "invalid level") ;
        }
    }

    //--------------------------------------------------------------------------
    // check the parent of each node
    //--------------------------------------------------------------------------

    if (parent_in != NULL)
    {
        for (int64_t i = 0 ; i < n ; i++)
        {
            if (i == src)
            {
                // src node is its own parent
                bool ok = (parent_in [src] == src) && (visited [src]) ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
            }
            else if (visited [i])
            {
                int64_t pi = parent_in [i] ;
                // ensure the parent pi is valid and has been visited
                bool ok = (pi >= 0 && pi < n) && visited [pi] ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
                // ensure the edge (pi,i) exists
                bool x ;
                int info = GrB_Matrix_extractElement_BOOL (&x, G->A, pi, i) ;
                ok = (info == GrB_SUCCESS) ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
                // ensure the parent's level is ok
                ok = (level_check [i] == level_check [pi] + 1) ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGraph_FREE_WORK ;

    if (print_timings)
    {
        LAGraph_Toc (&tt, tic, msg) ;
        printf ("LG_check_bfs check time: %g sec\n", tt) ;
    }
    return (0) ;
}

