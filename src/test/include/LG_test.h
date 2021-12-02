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
// LG_brutal_*:  brutal memory tests
//------------------------------------------------------------------------------

// Brutal memory tests use a global variable (LG_brutal) to tell them how many
// times malloc may be successfully called.  Once this counter reaches zero,
// LG_brutal_malloc, LG_brutal_calloc, and LG_brutal_realloc all return NULL.
// These tests must be used with care on methods that unpack/pack their input
// matrix G->A (such as several of the LG_check_* methods).  Those methods will
// return an empty G->A with no entries, if they fail in the middle.

// SuiteSparse:GraphBLAS is required for brutal memory testing.  The
// LG_brutal_malloc/calloc/realloc/free methods are passed to SuiteSparse
// GraphBLAS via GxB_init.  This way, out-of-memory conditions returned by
// GraphBLAS can be checked and handled by LAGraph.

// Use LG_brutal_setup to start LAGraph for brutal memory tests, and
// LG_brutal_teardown to finish.  To test a method with brutal memory tests,
// use LG_BRUTAL (LAGraph_method (...)).  The LAGraph_method(...) will be
// called with LG_brutal set to 0 (no mallocs allowed), 1 (one malloc), 2, ...
// until the method succeeds by returning a result >= 0.  If a method never
// returns a non-negative result, LG_BRUTAL will get stuck in an infinite loop.

// If LG_brutal starts as negative, then the brutal memory tests are not
// performed, and LG_brutal_malloc/calloc/etc never pretend to fail.

// A count (LG_nmalloc) is kept of the number of allocated blocks that have not
// yet been freed.  If this count is not zero after finalizing GraphBLAS
// and LAGraph, an error is reported.  No report is provided as to what blocks
// of memory are still allocated; use valgrind if that occurs.

// For methods that access files, or have other side effects, the LG_BRUTAL
// macro will not work, as LG_BRUTAL (LAGraph_MMRead (...)) for example.  Each
// iteration of the brutal loop must also close the file and reopen for the
// next brutal trial.

LAGRAPH_PUBLIC int LG_brutal_setup (char *msg) ;
LAGRAPH_PUBLIC int LG_brutal_teardown (char *msg) ;

LAGRAPH_PUBLIC int64_t LG_brutal ;
LAGRAPH_PUBLIC int64_t LG_nmalloc ;

LAGRAPH_PUBLIC
void *LG_brutal_malloc      // return pointer to allocated block of memory
(
    size_t size             // # of bytes to allocate
) ;

LAGRAPH_PUBLIC
void *LG_brutal_calloc      // return pointer to allocated block of memory
(
    size_t nitems,          // # of items to allocate
    size_t itemsize         // # of bytes per item
) ;

LAGRAPH_PUBLIC
void LG_brutal_free
(
    void *p                 // block to free
) ;

LAGRAPH_PUBLIC
void *LG_brutal_realloc     // return pointer to reallocated memory
(
    void *p,                // block to realloc
    size_t size             // new size of the block
) ;

// brutal memory testing of a GraphBLAS or LAGraph method, no burble
#define LG_BRUTAL(Method)                                       \
{                                                               \
    for (int nbrutal = 0 ; ; nbrutal++)                         \
    {                                                           \
        /* allow for only nbrutal mallocs before 'failing' */   \
        LG_brutal = nbrutal ;                                   \
        /* try the method with brutal malloc */                 \
        int brutal_result = Method ;                            \
        if (brutal_result >= 0)                                 \
        {                                                       \
            /* the method finally succeeded */                  \
            break ;                                             \
        }                                                       \
    }                                                           \
    LG_brutal = -1 ;  /* turn off brutal mallocs */             \
}

// brutal memory testing of a GraphBLAS or LAGraph method, and print results
#define LG_BRUTAL_BURBLE(Method)                                \
{                                                               \
    printf ("brutal test at line %4d: LG_nmalloc: %g\n",        \
        __LINE__, (double) LG_nmalloc) ;                        \
    printf ("method: " LG_XSTR(Method) "\n") ;                  \
    for (int nbrutal = 0 ; ; nbrutal++)                         \
    {                                                           \
        /* allow for only nbrutal mallocs before 'failing' */   \
        LG_brutal = nbrutal ;                                   \
        /* try the method with brutal malloc */                 \
        int brutal_result = Method ;                            \
        if (brutal_result >= 0)                                 \
        {                                                       \
            /* the method finally succeeded */                  \
            printf ("brutal test at line %4d: LG_nmalloc: %g,"  \
                " succeeded with %d mallocs\n", __LINE__,       \
                (double) LG_nmalloc, nbrutal) ;                 \
            break ;                                             \
        }                                                       \
    }                                                           \
    LG_brutal = -1 ;  /* turn off brutal mallocs */             \
}

#endif
