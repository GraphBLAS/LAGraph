//------------------------------------------------------------------------------
// LAGraph_tic:  a portable timer for accurate performance measurements
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

// There is no method that works on all operating systems for finding the
// current time with high resolution that is suitable for performance
// measurements.  The LAGraph_tic and LAGraph_toc functions provide a portable
// alternative.

// LAGraph_tic (tic) ; gets the current time and saves it in tic [0..1].

// t = LAGraph_toc (tic) ; returns the time in seconds since the last call to
// LAGraph_toc, as a single double value.

// Example:

/*
    double tic [2], t ;
    LAGraph_tic (tic) ;
    // ... do stuff
    t = LAGraph_toc (tic) ;
    printf ("time to 'do stuff' : %g (seconds)\n', t) ;
    // ... more stuff
    t = LAGraph_toc (tic) ;
    printf ("time to 'do stuff' and 'more stuff': %g (seconds)\n', t) ;
*/

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// LAGrapph_tic: return the current wallclock time in high resolution
//------------------------------------------------------------------------------

void LAGraph_tic            // returns current time in seconds and nanoseconds
(
    double tic [2]          // tic [0]: seconds, tic [1]: nanoseconds
)
{

    #if defined ( _OPENMP )

        // OpenMP is available; use the OpenMP timer function
        tic [0] = omp_get_wtime ( ) ;
        tic [1] = 0 ;

    #elif defined ( __linux__ )

        // Linux has a very low resolution clock() function, so use the high
        // resolution clock_gettime instead.  May require -lrt
        struct timespec t ;
        clock_gettime (CLOCK_MONOTONIC, &t) ;
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

}

