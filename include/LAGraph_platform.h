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
#include <GraphBLAS.h>

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

// what is the correct system 64-bit maximum integer?
#if defined(GxB_SUITESPARSE_GRAPHBLAS)
   #define LAGRAPH_INDEX_MAX GxB_INDEX_MAX
#else
   #define LAGRAPH_INDEX_MAX ULONG_MAX
#endif

#endif  // LAGRAPH_PLATFORM_H
