//------------------------------------------------------------------------------
// LAGraph_get_nthreads: get the # of threads to use
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// LAGraph_get_nthreads: get # of threads to use
// contributed by Tim Davis, Texas A&M

// See also LAGraph_get_nthreads

#include "LAGraph_internal.h"

int LAGraph_get_nthreads        // returns # threads to use, 1 if unknown
(
    void
)
{

    int nthreads = 1 ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    info = GxB_get (GxB_NTHREADS, &nthreads) ;
    if (info != GrB_SUCCESS)
    {
        fprintf (stderr, "LAGraph error:\n[%d]\n%s\nFile: %s Line: %d\n",
            info, GrB_error ( ), __FILE__, __LINE__) ;
        return (-9999) ;
    }
    #elif defined ( _OPENMP )
    nthreads = omp_get_max_threads ( ) ;
    #else
    // nothing to do ...
    #endif

    return (nthreads) ;
}

