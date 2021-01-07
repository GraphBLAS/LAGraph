//------------------------------------------------------------------------------
// LAGraph_Internal.h: include file for use within LAGraph itself
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

#ifndef LAGRAPH_INTERNAL_H
#define LAGRAPH_INTERNAL_H

//------------------------------------------------------------------------------
// include files
//------------------------------------------------------------------------------

#include "LAGraph2.h"

#if defined ( __linux__ )
#include <malloc.h>
#endif

//------------------------------------------------------------------------------
// LAGraph_CLEAR_MSG: clear the error msg string
//------------------------------------------------------------------------------

// When an LAGraph method starts, it first clears the caller's msg string.
#define LAGraph_CLEAR_MSG               \
{                                       \
    if (msg != NULL) msg [0] = '\0' ;   \
}

//------------------------------------------------------------------------------
// LAGraph_ERROR_MSG: set the error msg string
//------------------------------------------------------------------------------

// When an LAGraph method encounters an error, it can report details in the
// msg.  For example:

/*
    if (src < 0 || src >= n)
    {
        LAGraph_ERROR_MSG ("Source node %ld must be in range 0 to n-1,\n"
            "where n = %ld is the number of nodes in the graph.", src, n) ;
        return (-1) ;
    }
*/

#define LAGraph_ERROR_MSG(...)                                      \
{                                                                   \
    if (msg != NULL) snprintf (msg, LAGRAPH_MSG_LEN, __VA_ARGS__) ; \
}

//------------------------------------------------------------------------------
// LAGraph_FREE_WORK: free all workspace
//------------------------------------------------------------------------------

#ifndef LAGraph_FREE_WORK
#define LAGraph_FREE_WORK ;
#endif

//------------------------------------------------------------------------------
// LAGraph_FREE_ALL: free all workspace and all output arguments, on error
//------------------------------------------------------------------------------

#ifndef LAGraph_FREE_ALL
#define LAGraph_FREE_ALL        \
{                               \
    LAGraph_FREE_WORK ;         \
}
#endif

//------------------------------------------------------------------------------
// GrB_CATCH: catch an error from GraphBLAS
//------------------------------------------------------------------------------

// A simple GrB_CATCH macro to be used by GrB_TRY.  If the LAGraph function
// wants something else, then #define a GrB_CATCH macro before the #include
// "LAGraph_Internal.h" statement.

#ifndef GrB_CATCH
#define GrB_CATCH(info)                                 \
{                                                       \
    LAGraph_ERROR_MSG ("%s, line %d: failure: %d\n",    \
        __FILE__, __LINE__, info) ;                     \
    LAGraph_FREE_ALL ;                                  \
    return (-1) ;                                       \
}
#endif

//------------------------------------------------------------------------------
// LAGraph_CATCH: catch an error from LAGraph
//------------------------------------------------------------------------------

// A simple LAGraph_CATCH macro:  same as GrB_CATCH
#ifndef LAGraph_CATCH
#define LAGraph_CATCH(status) GrB_CATCH(status)
#endif

//------------------------------------------------------------------------------
// LAGraph_CHECK: check an error condition
//------------------------------------------------------------------------------

#define LAGraph_CHECK(error_condition,error_status,...) \
{                                                       \
    if (error_condition)                                \
    {                                                   \
        LAGraph_ERROR_MSG (__VA_ARGS__) ;               \
        LAGraph_FREE_ALL ;                              \
        return (error_status) ;                         \
    }                                                   \
}

//------------------------------------------------------------------------------
// LAGraph_CHECK_INIT: clear msg and do basic tests of a graph
//------------------------------------------------------------------------------

#define LAGraph_CHECK_INIT(G,msg)                                           \
{                                                                           \
    LAGraph_CLEAR_MSG ;                                                     \
    LAGraph_CHECK (G == NULL, -1, "graph is NULL") ;                        \
    LAGraph_CHECK (G->A == NULL, -1, "graph adjacency matrix is NULL") ;    \
    LAGraph_CHECK (G->kind <= LAGRAPH_UNKNOWN ||                            \
        G->kind > LAGRAPH_ADJACENCY_DIRECTED, -1, "graph kind invalid") ;   \
}

//------------------------------------------------------------------------------
// code development settings
//------------------------------------------------------------------------------

// turn off debugging; do not edit these three lines
#ifndef NDEBUG
#define NDEBUG
#endif

// These flags are used for code development.  Uncomment them as needed.

// to turn on debugging, uncomment this line:
// #undef NDEBUG

#undef ASSERT

#ifndef NDEBUG

    // debugging enabled
    #ifdef MATLAB_MEX_FILE
        // debugging when LAGraph is part of a mexFunction
        #define ASSERT(x)                                               \
        {                                                               \
            if (!(x)) mexErrMsgTxt ("failure: " __FILE__ " line: ") ;   \
        }
    #else
        #include <assert.h>
        #define ASSERT(x) assert (x) ;
    #endif

#else

    // debugging disabled
    #define ASSERT(x)

#endif

//------------------------------------------------------------------------------
// LAGraph_Multiply_size_t:  c = a*b but check for overflow
//------------------------------------------------------------------------------

static bool LAGraph_Multiply_size_t    // true if ok, false if overflow
(
    size_t *c,                  // c = a*b, or zero if overflow occurs
    const size_t a,
    const size_t b
)
{

    ASSERT (c != NULL) ;

    (*c) = 0 ;
    if (a == 0 || b == 0)
    { 
        return (true) ;
    }

    if (a > SIZE_MAX / 2 || b > SIZE_MAX / 2)
    { 
        // a or b are out of range
        return (false) ;
    }

    // a + b is now safe to compute
    if ((a + b) > (SIZE_MAX / LAGraph_MIN (a,b)))
    { 
        // a * b may overflow
        return (false) ;
    }

    // a * b will not overflow
    (*c) = a * b ;
    return (true) ;
}

//------------------------------------------------------------------------------
// Matrix Market format
//------------------------------------------------------------------------------

// %%MatrixMarket matrix <fmt> <type> <storage> uses the following enums:

typedef enum
{
    MM_coordinate,
    MM_array,
}
MM_fmt_enum ;

typedef enum
{
    MM_real,
    MM_integer,
    MM_complex,
    MM_pattern
}
MM_type_enum ;

typedef enum
{
    MM_general,
    MM_symmetric,
    MM_skew_symmetric,
    MM_hermitian
}
MM_storage_enum ;

// maximum length of each line in the Matrix Market file format

// The MatrixMarket format specificies a maximum line length of 1024.
// This is currently sufficient for GraphBLAS but will need to be relaxed
// if this function is extended to handle arbitrary user-defined types.
#define MMLEN 1024
#define MAXLINE MMLEN+6

#endif

