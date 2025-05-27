//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/kt_demo.c: test AllKTruss many times
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
    LAGraph_Free ((void **) &Cset, NULL) ;      \
    LAGraph_Free ((void **) &ntris, NULL) ;     \
    LAGraph_Free ((void **) &nedges, NULL) ;    \
    LAGraph_Free ((void **) &nsteps, NULL) ;    \
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

    GrB_Index n ;
    int64_t kmax ;
    GrB_Matrix_nrows (&n, G->A) ;
    LAGraph_Calloc ((void **) &Cset  , n, sizeof (GrB_Matrix), msg) ;
    LAGraph_Malloc ((void **) &ntris , n, sizeof (int64_t), msg) ;
    LAGraph_Malloc ((void **) &nedges, n, sizeof (int64_t), msg) ;
    LAGraph_Malloc ((void **) &nsteps, n, sizeof (int64_t), msg) ;

    LAGraph_AllKTruss (Cset, &kmax, ntris, nedges, nsteps, G, msg) ;
    int64_t kmax_ok = kmax ;

    //--------------------------------------------------------------------------
    // call AllKTruss many times
    //--------------------------------------------------------------------------

    #define NTRIALS 0
    printf ("AllKTruss: %d trials\n", NTRIALS) ;

    t = LAGraph_WallClockTime ( ) ;
    double t1 = t ;
    for (int k = 0 ; k < NTRIALS ; k++)
    {
        LAGraph_AllKTruss (Cset, &kmax, ntris, nedges, nsteps, G, msg) ;
        printf ("trial %d : all k-truss: kmax %g\n", k, (double) kmax) ;
        double tt = LAGraph_WallClockTime ( ) - t1 ;
        if (tt > 3)
        {
            printf ("%" PRId64 " : %d ok, %g sec\n", kmax, k,
                LAGraph_WallClockTime ( ) - t) ;
            fflush (stdout) ;
            t1 = LAGraph_WallClockTime ( ) ;
        }
        if (kmax != kmax_ok)
        {
            printf ("Abort! %" PRId64 " %" PRId64 "\n", kmax, kmax_ok) ;
            fflush (stdout) ;
            abort ( ) ;
        }
    }
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for %d trials:          %g sec\n", NTRIALS, t) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

