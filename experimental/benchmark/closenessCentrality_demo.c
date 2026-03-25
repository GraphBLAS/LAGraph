//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/closenessCentrality_demo.c: demo for Closeness
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Karan Bhalla and Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#define LG_FREE_ALL                         \
    {                                       \
        GrB_free (&centrality) ;            \
        GrB_free (&A) ;                     \
        GrB_free (&sources) ;               \
        LAGraph_Delete (&G, msg) ;          \
    }

#include "LAGraphX.h"
#include "LG_internal.h"
#include "LG_Xtest.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Vector centrality = NULL, sources = NULL ;

    LAGRAPH_TRY (LAGraph_Init (msg)) ;

    //--------------------------------------------------------------------------
    // parse command-line arguments
    //--------------------------------------------------------------------------

    // Usage:
    //   closenessCentrality_demo <matrix-market-file>
    //                            [num_sources] [use_weights] [use_floyd_warshall]
    //
    // <matrix-market-file>  : required; path to .mtx file
    // [num_sources]         : optional; number of random source nodes to score
    //                         (0 or omitted => score all nodes)
    // [use_weights]         : optional 0/1 (default 0); use edge weights
    // [use_floyd_warshall]  : optional 0/1 (default 0); use Floyd-Warshall APSP
    //                         (only valid when num_sources is 0/omitted)

    if (argc < 2 || argc > 5)
    {
        printf ("Usage: %s <matrix-market-file>"
                " [num_sources] [use_weights] [use_floyd_warshall]\n",
                argv [0]) ;
        return (GrB_INVALID_VALUE) ;
    }

    bool use_weights        = false ;
    bool use_floyd_warshall = false ;
    int  num_sources        = 0 ;

    if (argc > 2) num_sources        = atoi (argv [2]) ;
    if (argc > 3) use_weights        = (atoi (argv [3]) != 0) ;
    if (argc > 4) use_floyd_warshall = (atoi (argv [4]) != 0) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    FILE *f = fopen (argv [1], "r") ;
    if (f == NULL)
    {
        printf ("Error: unable to open file %s\n", argv [1]) ;
        LG_FREE_ALL ;
        return (GrB_INVALID_VALUE) ;
    }

    double t = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGraph_MMRead (&A, f, msg)) ;
    fclose (f) ;

    uint64_t n ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    LAGRAPH_TRY (LAGraph_DeleteSelfEdges (G, msg)) ;
    LAGRAPH_TRY (LAGraph_Cached_AT (G, msg)) ;

    if (use_weights)
    {
        LAGRAPH_TRY (LAGraph_Cached_EMin (G, msg)) ;
    }

    t = LAGraph_WallClockTime () - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // build optional source-node vector
    //--------------------------------------------------------------------------

    if (num_sources > 0)
    {
        if (use_floyd_warshall)
        {
            printf ("Warning: use_floyd_warshall ignored when num_sources > 0;"
                    " falling back to per-node shortest paths.\n") ;
            use_floyd_warshall = false ;
        }

        if ((uint64_t) num_sources >= n)
        {
            printf ("Error: num_sources (%" PRId32 ") must be less than"
                    " n (%" PRIu64 ").\n", num_sources, n) ;
            LG_FREE_ALL ;
            return (GrB_INVALID_VALUE) ;
        }

        GRB_TRY (GrB_Vector_new (&sources, GrB_UINT64, num_sources)) ;

        double t_seed = LAGraph_WallClockTime () ;
        srand ((int) t_seed) ;

        bool *used = NULL ;
        LAGRAPH_TRY (LAGraph_Calloc ((void **) &used, n, sizeof (bool), msg)) ;

        for (int i = 0 ; i < num_sources ; i++)
        {
            GrB_Index idx ;
            do { idx = rand () % n ; } while (used [idx]) ;
            used [idx] = true ;
            GRB_TRY (GrB_Vector_setElement (sources, idx, i)) ;
        }
        LAGRAPH_TRY (LAGraph_Free ((void **) &used, msg)) ;

        printf ("Using %d random source nodes.\n", num_sources) ;
    }

    //--------------------------------------------------------------------------
    // compute closeness centrality
    //--------------------------------------------------------------------------

    printf ("\nCloseness params: use_weights=%d use_floyd_warshall=%d"
            " num_sources=%d\n\n",
            (int) use_weights, (int) use_floyd_warshall, num_sources) ;

    t = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGr_ClosenessCentrality (&centrality, G, sources,
                                           use_weights, use_floyd_warshall,
                                           msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("Time for LAGr_ClosenessCentrality: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // print results
    //--------------------------------------------------------------------------

    GrB_Index nvals = 0 ;
    GRB_TRY (GrB_Vector_nvals (&nvals, centrality)) ;
    printf ("Number of scored nodes: %" PRIu64 "\n", (uint64_t) nvals) ;

    // print first few scores
    GrB_Index print_limit = (n < 20) ? n : 20 ;
    printf ("Closeness centrality (first %" PRIu64 " nodes):\n",
            (uint64_t) print_limit) ;
    for (GrB_Index i = 0 ; i < print_limit ; i++)
    {
        double score = 0 ;
        GrB_Info info = GrB_Vector_extractElement_FP64 (&score, centrality, i) ;
        if (info == GrB_SUCCESS)
        {
            printf ("  node %" PRIu64 ": %g\n", (uint64_t) i, score) ;
        }
        else
        {
            printf ("  node %" PRIu64 ": (unreachable)\n", (uint64_t) i) ;
        }
    }

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
