//------------------------------------------------------------------------------
// LAGraph_log:  log a test result for LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite

// LAGraph_log:  log a result from LAGraph, contributed by Tim Davis, Texas A&M

#include <time.h>
#include <LAGraph.h>
#include <LAGraphX.h>

#define LAGraph_FREE_ALL

//****************************************************************************
static void delete_newline (char *p)
{
    for ( ; *p ; p++)
    {
        if (*p == '\n')
        {
            *p = '\0' ;
            return ;
        }
    }
}

//****************************************************************************
GrB_Info LAGraph_log
(
    char *caller,           // calling function
    char *message1,         // message to include (may be NULL)
    char *message2,         // message to include (may be NULL)
    int nthreads,           // # of threads used
    double t                // time taken by the test
)
{

    GrB_Info info ;

    // try to get the hostname and cpu type
    char hostname [1024] ;
    hostname [0] = '\0' ;
    FILE *host = fopen ("/etc/hostname", "r") ;
    if (host != NULL)
    {
        // get the hostname
        if (!fgets (hostname, 1000, host)) hostname [0] = '\0' ;
        // delete the newline
        delete_newline (hostname) ;
        fclose (host) ;
    }

    char filename [2048] ;
    snprintf (filename, 2000, "log_%s.txt", hostname) ;

    FILE *f = fopen (filename, "a") ;

    if (f == NULL)
    {
        LAGRAPH_ERROR ("cannot open logfile", GrB_INVALID_VALUE) ;
    }

    time_t now = time (NULL) ;
    fprintf (f, "\nFrom: %s\nDate: %s", caller, ctime (&now)) ;

    FILE *cpuinfo = fopen ("/proc/cpuinfo", "r") ;
    if (cpuinfo != NULL)
    {
        while (1)
        {
            // get the line
            char line [2048] ;
            if (!fgets (line, 2000, cpuinfo))
            {
                break ;
            }
            // find the colon
            bool found = false ;
            char *p = line ;
            for ( ; (*p) ; p++)
            {
                if (*p == ':')
                {
                    *p = '\0' ;
                    p++ ;
                    if (strncmp (line, "model name", 10) == 0)
                    {
                        delete_newline (p) ;
                        fprintf (f, "CPU: %s\n", p) ;
                        found = true ;
                        break ;
                    }
                }
            }
            if (found) break ;
        }
    }

    int max_threads;
    LAGraph_GetNumThreads (&max_threads, NULL);
    fprintf (f, "max # of threads: %d\n", max_threads) ;

#if LG_SUITESPARSE
    char *library_date ;
    GxB_get (GxB_LIBRARY_DATE, &library_date) ;
    fprintf (f, "SuiteSparse:GraphBLAS %s\n", library_date) ;
#endif

    fprintf (f, "Message: %s : %s\n# threads used: %d time: %g\n",
        (message1 == NULL) ? "" : message1,
        (message2 == NULL) ? "" : message2,
        nthreads, t) ;

    fclose (f) ;
    return (GrB_SUCCESS) ;
}
