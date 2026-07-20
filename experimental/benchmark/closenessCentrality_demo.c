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

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include <stdio.h>
#include <stdlib.h>

#define LG_FREE_ALL                         \
    {                                       \
        GrB_free (&centrality) ;            \
        GrB_free (&sources) ;               \
        GrB_free (&Delta) ;                 \
        LAGraph_Delete (&G, msg) ;          \
    }

int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    LAGraph_Graph G = NULL ;
    GrB_Vector centrality = NULL, sources = NULL ;
    GrB_Scalar Delta = NULL ;

    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // parse command-line arguments
    //--------------------------------------------------------------------------

    // Usage:
    //   closenessCentrality_demo <matrix-market-file> <algorithm>
    //                            [num_sources] [use_weights] [delta]
    //
    // <matrix-market-file>  : required; path to .mtx file
    // <algorithm>           : required int; shortest-path algorithm:
    //                           0 = CC_BFS           (unweighted BFS)
    //                           1 = CC_SSSP          (delta-stepping SSSP)
    //                           2 = CC_BELLMAN_FORD  (Bellman-Ford)
    //                           3 = CC_FLOYD_WARSHALL (all-pairs, no sources)
    // [num_sources]         : optional; number of random source nodes to score
    //                         (0 or omitted => score all nodes)
    // [use_weights]         : optional 0/1 (default 0); use edge weights
    // [delta]               : optional float; delta for CC_SSSP
    //                         (ignored for other algorithms)

    if (argc < 3 || argc > 6)
    {
        printf ("Usage: %s <matrix-market-file> <algorithm>"
                " [num_sources] [use_weights] [delta]\n"
                "  algorithm: 0=CC_BFS  1=CC_SSSP  2=CC_BELLMAN_FORD"
                "  3=CC_FLOYD_WARSHALL\n",
                argv [0]) ;
        return (GrB_INVALID_VALUE) ;
    }

    LAGraph_cc_algo_t  algorithm   = (LAGraph_cc_algo_t) atoi (argv [2]) ;
    bool       use_weights = false ;
    double     delta_val   = -1 ;
    int        num_sources = 0 ;

    if (argc > 3) num_sources = atoi (argv [3]) ;
    if (argc > 4) use_weights = (atoi (argv [4]) != 0) ;
    if (argc > 5) delta_val   = strtod (argv [5], NULL) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    // char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    double t = LAGraph_WallClockTime () ;
    LAGRAPH_TRY(
        readproblem(&G, NULL, false, false, false, NULL, false, argc, argv)) ;

    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    LAGRAPH_TRY (LAGraph_Cached_AT (G, msg)) ;

    if (use_weights)
    {
        LAGRAPH_TRY (LAGraph_Cached_EMin (G, msg)) ;
    }

    GrB_Index nvals ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    t = LAGraph_WallClockTime () - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;
    printf ("Nodes: %" PRIu64 "  Edges: %" PRIu64 "\n",
            (uint64_t) n, (uint64_t) nvals) ;

    printf ("\n==========================\nThe input graph matrix G:\n") ;
    // LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // build optional source-node vector
    //--------------------------------------------------------------------------

    if (num_sources > 0)
    {
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

    printf ("\nCloseness params: use_weights=%d algorithm=%d"
            " num_sources=%d\n\n",
            (int) use_weights, (int) algorithm, num_sources) ;

    if (algorithm == CC_SSSP && delta_val > 0)
    {
        GRB_TRY (GrB_Scalar_new (&Delta, GrB_FP64)) ;
        GRB_TRY (GrB_Scalar_setElement_FP64 (Delta, delta_val)) ;
    }

    t = LAGraph_WallClockTime () ;
    LAGRAPH_TRY (LAGr_ClosenessCentrality (&centrality, G, sources,
                                           use_weights, algorithm, Delta,
                                           msg)) ;
    t = LAGraph_WallClockTime () - t ;
    printf ("Time for LAGr_ClosenessCentrality: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // print results
    //--------------------------------------------------------------------------

    GrB_Index scoredvals = 0 ;
    GRB_TRY (GrB_Vector_nvals (&scoredvals, centrality)) ;
    printf ("Number of scored nodes: %" PRIu64 "\n", (uint64_t) scoredvals) ;

    // print first few scores
    // GrB_Index print_limit = (n < 20) ? n : 20 ;
    // printf ("Closeness centrality (first %" PRIu64 " nodes):\n",
    //         (uint64_t) print_limit) ;
    // for (GrB_Index i = 0 ; i < print_limit ; i++)
    // {
    //     double score = 0 ;
    //     GrB_Info info = GrB_Vector_extractElement_FP64 (&score, centrality, i) ;
    //     if (info == GrB_SUCCESS)
    //     {
    //         printf ("  node %" PRIu64 ": %g\n", (uint64_t) i, score) ;
    //     }
    //     else
    //     {
    //         printf ("  node %" PRIu64 ": (unreachable)\n", (uint64_t) i) ;
    //     }
    // }

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
