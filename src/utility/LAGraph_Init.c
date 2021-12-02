//------------------------------------------------------------------------------
// LAGraph_Init: start GraphBLAS and LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Init (char *msg)        // return 0 if success, -1 if failure
{

    LG_CLEAR_MSG ;

    #ifdef __linux__
    // Use mallopt to speedup malloc and free on Linux.  Otherwise, it can take
    // several seconds to free a large block of memory.  For this to be
    // effective, LAGraph_Init must be called before the user program does any
    // mallocs or frees itself.
    mallopt (M_MMAP_MAX, 0) ;           // disable mmap; it's too slow
    mallopt (M_TRIM_THRESHOLD, -1) ;    // disable sbrk trimming
    mallopt (M_TOP_PAD, 16*1024*1024) ; // increase padding to speedup malloc
    #endif

    // use ANSI C memory allocation functions
    int result = LAGraph_Xinit (malloc, calloc, realloc, free, msg) ;
    LG_CHECK (result != GrB_SUCCESS, result, "failed to initialize LAGraph") ;
    return (result) ;
}
