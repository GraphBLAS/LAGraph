//------------------------------------------------------------------------------
// LAGraph_Tic: start the timer
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// Example usage:

/*
    double tic [2], t ;
    LAGraph_Tic (tic, NULL) ;
    // ... do stuff
    LAGraph_Toc (&t, tic, NULL) ;
    printf ("time to 'do stuff' : %g (seconds)\n', t) ;
    // ... more stuff
    LAGraph_Toc (&t, tic, NULL) ;
    printf ("time to 'do stuff' and 'more stuff': %g (seconds)\n', t) ;
*/

#include "LG_internal.h"

#if !defined ( _OPENMP )
    #include <time.h>
    #if defined ( __linux__ ) || defined ( __GNU__ )
        #include <sys/time.h>
    #endif
    #if defined ( __MACH__ ) && defined ( __APPLE__ )
        #include <mach/clock.h>
        #include <mach/mach.h>
    #endif
#endif

int LAGraph_Tic
(
    double tic [2],     // tic [0]: seconds, tic [1]: nanoseconds
    char *msg
)
{

    LG_CLEAR_MSG ;

    #if defined ( _OPENMP )

        // OpenMP is available; use the OpenMP timer function
        tic [0] = omp_get_wtime ( ) ;
        tic [1] = 0 ;

    #elif defined ( __linux__ )

        // Linux has a very low resolution clock() function, so use the high
        // resolution clock_gettime instead.  May require -lrt
        struct timespec t ;
        LAGRAPH_TRY (clock_gettime (CLOCK_MONOTONIC, &t)) ;
        tic [0] = (double) t.tv_sec ;
        tic [1] = (double) t.tv_nsec ;

    #elif defined ( __MACH__ )

        // Mac OSX
        clock_serv_t cclock ;
        mach_timespec_t t ;
        host_get_clock_service (mach_host_self ( ), SYSTEM_CLOCK, &cclock) ;
        clock_get_time (cclock, &t) ;
        mach_port_deallocate (mach_task_self ( ), cclock) ;
        tic [0] = (double) t.tv_sec;
        tic [1] = (double) t.tv_nsec;

    #else

        // The ANSI C11 clock() function is used instead.  This gives the
        // processor time, not the wallclock time, and it might have low
        // resolution.  It returns the time since some unspecified fixed time
        // in the past, as a clock_t integer.  The clock ticks per second are
        // given by CLOCKS_PER_SEC.  In Mac OSX this is a very high resolution
        // clock, and clock ( ) is faster than clock_get_time (...) ;
        clock_t t = clock ( ) ;
        tic [0] = ((double) t) / ((double) CLOCKS_PER_SEC) ;
        tic [1] = 0 ;

    #endif

    return (GrB_SUCCESS) ;
}

