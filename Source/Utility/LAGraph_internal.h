//------------------------------------------------------------------------------
// LAGraph_internal.h: include file for internal use in LAGraph
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// Include file for LAGraph functions.  This file should not be included in
// user applications.  See LAGraph.h instead.

#include "LAGraph.h"
#include <ctype.h>

#ifndef CMPLX
#define CMPLX(real,imag) \
    ( \
    (double complex)((double)(real)) + \
    (double complex)((double)(imag) * _Complex_I) \
    )
#endif

// "I" is used in LAGraph and GraphBLAS to denote a list of row indices; remove
// it here
#undef I

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
// LAGRAPH_OK: call LAGraph or GraphBLAS and check the result
//------------------------------------------------------------------------------

// The including file must define LAGRAPH_FREE_ALL as a macro that frees all
// workspace if an error occurs.  method can be a GrB_Info scalar as well,
// so that LAGRAPH_OK(info) works.

#define LAGRAPH_OK(method)                                                  \
{                                                                           \
    GrB_Info this_info = method ;                                           \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))          \
    {                                                                       \
        fprintf (stderr, "LAGraph error: [%d]\n%s\n",                       \
            this_info, GrB_error ( )) ;                                     \
        LAGRAPH_FREE_ALL ;                                                  \
        return (this_info) ;                                                \
    }                                                                       \
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

