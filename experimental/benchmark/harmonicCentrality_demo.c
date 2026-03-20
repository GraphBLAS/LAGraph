//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/harmonicCentrality_demo.c: benchmark for
// harmonic centrality
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel A. Gomez, FalkorDB

//------------------------------------------------------------------------------

// Usage:  harmonicCentrality_demo matrixmarketfile.mtx
//         harmonicCentrality_demo matrixmarketfile.mtx exact

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

#define LG_FREE_ALL                         \
{                                           \
    GrB_free (&A) ;                         \
    GrB_free (&scores_approx) ;             \
    GrB_free (&scores_exact) ;              \
    GrB_free (&node_weights) ;              \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;

    GrB_Matrix A = NULL ;
    LAGraph_Graph G = NULL ;
    GrB_Vector scores_approx = NULL ;
    GrB_Vector scores_exact = NULL ;
    GrB_Vector node_weights = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    if (argc < 2)
    {
        printf ("Usage: %s <matrix-market-file> [exact]\n", argv [0]) ;
        return (GrB_INVALID_VALUE) ;
    }

    bool run_exact = (argc >= 3 && strcmp (argv [2], "exact") == 0) ;

    char *matrix_name = argv [1] ;
    FILE *f = fopen (matrix_name, "r") ;
    if (f == NULL)
    {
        printf ("Error: unable to open file %s\n", matrix_name) ;
        return (GrB_INVALID_VALUE) ;
    }

    double t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_MMRead (&A, f, msg)) ;
    fclose (f) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    printf ("Graph: %s (%" PRIu64 " nodes, %" PRIu64 " edges)\n",
        matrix_name, (uint64_t) n, (uint64_t) nvals) ;

    // construct a graph
    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;

    //--------------------------------------------------------------------------
    // create boolean node_weights (all nodes participate with weight = 1)
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&node_weights, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_assign_BOOL (
        node_weights, NULL, NULL, true, GrB_ALL, n, NULL)) ;

    //--------------------------------------------------------------------------
    // compute approximate harmonic centrality (HLL)
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGr_HarmonicCentrality (
        &scores_approx, NULL, G, node_weights, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGr_HarmonicCentrality (approx): %g sec\n", t) ;

    // print a summary of the scores
    LAGraph_PrintLevel pr = (n <= 40) ? LAGraph_COMPLETE : LAGraph_SHORT ;
    printf ("\napproximate scores:\n") ;
    LAGRAPH_TRY (LAGraph_Vector_Print (scores_approx, pr, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // optionally compute exact harmonic centrality (BFS)
    //--------------------------------------------------------------------------

    if (run_exact)
    {
        t = LAGraph_WallClockTime ( ) ;
        LAGRAPH_TRY (LAGr_HarmonicCentrality_exact (
            &scores_exact, NULL, G, node_weights, node_weights, msg)) ;
        t = LAGraph_WallClockTime ( ) - t ;
        printf ("\nTime for LAGr_HarmonicCentrality_exact:  %g sec\n", t) ;

        printf ("\nexact scores:\n") ;
        LAGRAPH_TRY (LAGraph_Vector_Print (scores_exact, pr, stdout, msg)) ;
    }

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGraph_Delete (&G, msg) ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
