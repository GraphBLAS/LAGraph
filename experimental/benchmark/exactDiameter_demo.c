//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/exactDiameter_demo.c: a simple demo
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// This main program is a simple driver for testing and benchmarking the
// LAGraph_ExactDiameter algorithm, in experimental/algorithm based on 
// helloworld2_demo.  To use it, compile LAGraph while in the build folder
// with these commands:
//
//      cd LAGraph/build
//      cmake ..
//      make -j8
//
// Then run this demo with an input matrix.  For example:
//
//      ./experimental/benchmark/exactDiameter_demo < ../data/west0067.mtx
//      ./experimental/benchmark/exactDiameter_demo < ../data/karate.mtx
//      set up for raid graphs:
//      ./experimental/benchmark/exactDiameter_demo < /raid/matrices/com-Youtube/com-Youtube.mtx
//
#include "LAGraphX.h"

// LAGRAPH_CATCH is required by LAGRAPH_TRY.  If an error occurs, this macro
// catches it and takes corrective action, then terminates this program.
#define LAGRAPH_CATCH(info)                     \
{                                               \
    GrB_free (&Y) ;                             \
    GrB_free (&A) ;                             \
    GrB_free (&peripheral) ;                         \
    LAGraph_Delete (&G, msg) ;                  \
    return (info) ;                             \
}

// GRB_CATCH is required by GRB_TRY (although GRB_TRY isn't used here)
#define GRB_CATCH(info) LAGRAPH_CATCH(info)

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix Y = NULL, A = NULL ;
    GrB_Index diameter;
    GrB_Vector peripheral = NULL ;
    int numInBatch = 10;

    // start GraphBLAS and LAGraph
    LAGRAPH_TRY (LAGraph_Init (msg)) ;

    //--------------------------------------------------------------------------
    // read in the graph via a Matrix Market file from stdin
    //--------------------------------------------------------------------------

    double t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_MMRead (&A, stdin, msg)) ;
    LAGRAPH_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to read the graph:      %g sec\n", t) ;

    printf ("\n==========================The input graph matrix G:\n") ;
    LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // try the ExactDiameter algorithm
    //--------------------------------------------------------------------------

    printf ("\n==========================Running diameter\n") ;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_ExactDiameter (&diameter, &peripheral, NULL, G, numInBatch, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_ExactDiameter: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results (make sure outputs match)
    //--------------------------------------------------------------------------

    

    //--------------------------------------------------------------------------
    // print the results 
    //--------------------------------------------------------------------------

    int d = diameter;
    printf ("\n===============================Diameter found: %d \n", d) ;
    printf ("\n===============================The result peripheral vector is:\n") ;
    LAGRAPH_TRY (LAGraph_Vector_Print (peripheral, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    GrB_free (&Y) ;
    LAGraph_Delete (&G, msg) ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

