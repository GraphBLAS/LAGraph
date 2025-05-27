//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/rcc_demo.c: a demo for RichClubCoefficient
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

//      ./experimental/benchmark/rcc_demo ../data/west0067.mtx
//      ./experimental/benchmark/rcc_demo < ../data/west0067.mtx
//      ./experimental/benchmark/rcc_demo ../data/karate.mtx
//

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"

// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    LAGraph_Delete (&G, msg) ;                  \
    LAGraph_Delete (&G_new, msg) ;              \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL, G_new = NULL;
    double pQ = 0;

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
    GrB_Index swaps = (argc > 2) ? atoi(argv [2]): 100;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    FILE *f = (argc > 3) ? fopen (argv[3], "w"): NULL;
    LG_TRY (readproblem (
        &G,         // the graph that is read from stdin or a file
        NULL,       // source nodes (none, if NULL)
        true,      // make the graph undirected, if true
        true,      // remove self-edges, if true
        true,      // return G->A as structural, if true,
        GrB_BOOL,       // prefered GrB_Type of G->A; null if no preference
        false,      // ensure all entries are positive, if true
        argc, argv)) ;  // input to this main program
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LG_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // try the LAGraph_EdgeSwap algorithm
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Cached_OutDegree (G, msg)) ;
    printf("Time To Swap #################################################") ;
    t = LAGraph_WallClockTime ( ) ;
    LG_TRY (LAGraph_SwapEdges (&G_new, &pQ, G, swaps, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("LAGraph_SwapEdges took %g sec and made %g swaps\n", t, pQ) ;
    
    //--------------------------------------------------------------------------
    // check the results 
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    bool result = false;
    LG_TRY (LAGraph_Cached_OutDegree (G_new, msg)) ;
    LG_TRY (LAGraph_Vector_IsEqual(
        &result, G->out_degree, G_new->out_degree, msg)) ;
    if (result)
    {
        printf ("Test passed.\n") ;
    }
    else
    {
        printf ("Test failure!\n") ;
    }

    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf ("\n===============================The result matrix:\n") ;
    LG_TRY (LAGraph_Matrix_Print (G_new -> A, LAGraph_SHORT, stdout, msg)) ;
    if(f)
    {
        LG_TRY (binwrite(&(G_new -> A), f, NULL)) ;
    }
    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
