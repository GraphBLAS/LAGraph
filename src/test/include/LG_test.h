//------------------------------------------------------------------------------
// LG_test.h: include file for LAGraph test library
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#ifndef LG_TEST_H
#define LG_TEST_H

#include <GraphBLAS.h>
#include <LAGraph.h>

int LG_check_bfs
(
    // input
    GrB_Vector Level,       // optional; may be NULL
    GrB_Vector Parent,      // optional; may be NULL
    LAGraph_Graph G,
    GrB_Index src,
    char *msg
) ;

int LG_check_tri        // -1 if out of memory, 0 if successful
(
    // output
    uint64_t *ntri,     // # of triangles in A
    // input
    LAGraph_Graph G,    // the structure of G->A must be symmetric
    char *msg
) ;

int LG_check_cc
(
    // input
    GrB_Vector Component,   // Component(i)=k if node is in the kth Component
    LAGraph_Graph G,
    char *msg
) ;

bool LG_get_vector
(
    int64_t *x,
    GrB_Vector X,
    int64_t n,
    int64_t missing
) ;

int LG_check_sssp
(
    // input
    GrB_Vector Path_Length,     // Path_Length(i) is the length of the
                                // shortest path from src to node i.
    LAGraph_Graph G,            // all edge weights must be > 0
    GrB_Index src,
    char *msg
) ;

int LG_check_export
(
    // input
    LAGraph_Graph G,        // export G->A in CSR format
    // output
    GrB_Index **Ap_handle,  // size Ap_len on output
    GrB_Index **Aj_handle,  // size Aj_len on output
    void **Ax_handle,       // size Ax_len * typesize on output
    GrB_Index *Ap_len,
    GrB_Index *Aj_len,
    GrB_Index *Ax_len,
    size_t *typesize,       // size of the type of A
    char *msg
) ;

//------------------------------------------------------------------------------
// LG_check_malloc.c:  brutal memory tests
//------------------------------------------------------------------------------

LAGRAPH_PUBLIC int64_t LG_brutal ;
LAGRAPH_PUBLIC int64_t LG_nmalloc ;

LAGRAPH_PUBLIC
void *LG_check_malloc       // return pointer to allocated block of memory
(
    size_t size             // # of bytes to allocate
) ;

LAGRAPH_PUBLIC
void *LG_check_calloc       // return pointer to allocated block of memory
(
    size_t nitems,          // # of items to allocate
    size_t itemsize         // # of bytes per item
) ;

LAGRAPH_PUBLIC
void LG_check_free
(
    void *p                 // block to free
) ;

LAGRAPH_PUBLIC
void *LG_check_realloc      // return pointer to reallocated memory
(
    void *p,                // block to realloc
    size_t size             // new size of the block
) ;

#endif
