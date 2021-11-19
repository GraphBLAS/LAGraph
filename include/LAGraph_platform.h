//------------------------------------------------------------------------------
// LAGraph_platform.h: wrappers around platform-specific GraphBLAS extensions
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Abstract away as much about implementation specific aspects of GraphBLAS
// distributions and operating systems as possible.

//------------------------------------------------------------------------------

#ifndef LAGRAPH_PLATFORM_H
#define LAGRAPH_PLATFORM_H

//==============================================================================
// include files
//==============================================================================

#include <limits.h>

//==============================================================================

#if defined ( _OPENMP )
    #include <omp.h>
#endif

#if defined ( __linux__ ) || defined ( __GNU__ )
    #include <sys/time.h>
#endif

#if defined ( __MACH__ ) && defined ( __APPLE__ )
    #include <mach/clock.h>
    #include <mach/mach.h>
#endif

#if !defined(__cplusplus)
    #define LAGRAPH_RESTRICT restrict
#elif defined(_MSC_BUILD) || defined(__clang__) || defined(__GNUC__) || defined(__INTEL_COMPILER)
    #define LAGRAPH_RESTRICT __restrict
#else
    #define LAGRAPH_RESTRICT
#endif

//==============================================================================
// GraphBLAS platform specifics

// vanilla vs SuiteSparse:
#if !defined ( LG_VANILLA )
// by default, set LG_VANILLA to false
#define LG_VANILLA 0
#endif

#if ( !LG_VANILLA ) && defined ( GxB_SUITESPARSE_GRAPHBLAS )
    // use SuiteSparse, and its GxB* extensions
    #define LG_SUITESPARSE 1
#else
    // use any GraphBLAS library (possibly SuiteSparse) but with no GxB*
    #define LG_SUITESPARSE 0
#endif

#endif  // LAGRAPH_PLATFORM_H
