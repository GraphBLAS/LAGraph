//------------------------------------------------------------------------------
// LAGraph_xinit:  start GraphBLAS and LAGraph
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Initialize GraphBLAS, and then LAGraph.

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL                    \
{                                           \
    LAGraph_free_global ( ) ;               \
}

/*
GrB_Info GxB_init           // start up GraphBLAS and also define malloc, etc
(
    const GrB_Mode mode,    // blocking or non-blocking mode

    // pointers to memory management functions
    void * (* user_malloc_function  ) (size_t),
    void * (* user_calloc_function  ) (size_t, size_t),
    void * (* user_realloc_function ) (void *, size_t),
    void   (* user_free_function    ) (void *),

    // This argument appears in SuiteSparse:GraphBLAS 3.0.0 but not 2.x:
    bool user_malloc_is_thread_safe
) ;
*/

GrB_Info LAGraph_xinit
(
    // pointers to memory management functions
    void * (* user_malloc_function  ) (size_t),
    void * (* user_calloc_function  ) (size_t, size_t),
    void * (* user_realloc_function ) (void *, size_t),
    void   (* user_free_function    ) (void *),

    bool user_malloc_is_thread_safe
)
{

    // initialize GraphBLAS
    GrB_Info info ;

#if defined ( GxB_SUITESPARSE_GRAPHBLAS )

    #if ( GxB_IMPLEMENTATION_MAJOR >= 3 )

    LAGRAPH_OK (GxB_init (GrB_NONBLOCKING,
        user_malloc_function,
        user_calloc_function,
        user_realloc_function,
        user_free_function,
        user_malloc_is_thread_safe)) ;

    #else

    LAGRAPH_OK (GxB_init (GrB_NONBLOCKING,
        user_malloc_function,
        user_calloc_function,
        user_realloc_function,
        user_free_function)) ;

    #endif

    // save the memory management pointers in global LAGraph space
    LAGraph_malloc_function  = user_malloc_function ;
    LAGraph_calloc_function  = user_calloc_function ;
    LAGraph_realloc_function = user_realloc_function ;
    LAGraph_free_function    = user_free_function ;
    LAGraph_malloc_is_thread_safe = user_malloc_is_thread_safe ;

    // allocate all global objects
    LAGRAPH_OK (LAGraph_alloc_global ( )) ;

    return (GrB_SUCCESS) ;

#else

    // GxB_init is not available.  Use LAGraph_init instead.
    return (GrB_PANIC) ;

#endif
}

