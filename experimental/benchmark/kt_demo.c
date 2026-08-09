//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/kt_demo.c: test KTruss many times
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_Matrix_free (&C) ;                      \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix C = NULL, R = NULL ;
    GrB_Matrix *Cset = NULL ;
    int64_t *ntris = NULL, *nedges = NULL, *nsteps = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;               // set true for diagnostic outputs
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).

    double t = LAGraph_WallClockTime ( ) ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        true,       // make the graph undirected, if true
        true,       // remove self-edges, if true
        false,      // return G->A as structural, if true,
        NULL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

//  printf ("\n==========================The input graph matrix G:\n") ;
//  LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // compute each k-truss
    //--------------------------------------------------------------------------

    LAGraph_KTruss (&C, G, 3, msg) ;
    int64_t ntriangles = 0 ;
    GrB_Matrix_reduce_INT64 (&ntriangles, NULL, GrB_PLUS_MONOID_INT64, C, NULL) ;
    ntriangles = ntriangles / 6 ;
    printf ("# triangles: %ld\n", ntriangles) ;

    GrB_Matrix_free (&C) ;

    //--------------------------------------------------------------------------
    // call KTruss many times
    //--------------------------------------------------------------------------

    #define NTRIALS 1
    printf ("KTruss: %d trials\n", NTRIALS) ;

    t = LAGraph_WallClockTime ( ) ;
    for (int k = 0 ; k < NTRIALS ; k++)
    {
        GrB_Matrix_free (&C) ;
        LAGraph_KTruss (&C, G, 3, msg) ;
    }
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for %d trials:          %g sec\n", NTRIALS, t) ;

    uint64_t nvals ;
    GrB_Matrix_nvals (&nvals, C) ;
    printf ("Avg time: %g sec, 3-truss nvals: %lu\n", t / NTRIALS, nvals) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

