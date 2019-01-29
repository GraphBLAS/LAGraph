//------------------------------------------------------------------------------
// LAGraph_internal.h: include file for internal use in LAGraph
//------------------------------------------------------------------------------

// LAGraph, (TODO list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// This file should not be included in user applications.  See LAGraph.h
// instead.

#include "LAGraph.h"
#include <ctype.h>

#ifdef MATLAB_MEX_FILE
// compiling LAGraph for a MATLAB mexFunction.  Use mxMalloc, mxFree, etc.
#include "mex.h"
#include "matrix.h"
#define malloc  mxMalloc
#define free    mxFree
#define calloc  mxCalloc
#define realloc mxRealloc
#endif

//------------------------------------------------------------------------------
// code development settings
//------------------------------------------------------------------------------

// LAGRAPH_XSTR: convert the content of x into a string "x"
#define LAGRAPH_XSTR(x) LAGRAPH_STR(x)
#define LAGRAPH_STR(x) #x

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
    #define ASSERT(x) \
    {                                                                       \
        if (!(x))                                                           \
        {                                                                   \
            mexErrMsgTxt ("failure: " __FILE__ " line: "                    \
                LAGRAPH_XSTR(__LINE__)) ;                                   \
        }                                                                   \
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

