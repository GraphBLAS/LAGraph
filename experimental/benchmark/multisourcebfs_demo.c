//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/multisourcebfs_demo.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// This main program is a simple driver for testing and benchmarking the
// LAGraph_MultiSourceBFS algorithm, in experimental/algorithm based on 
// helloworld2_demo.  To use it, compile LAGraph while in the build folder
// with these commands:
//
//      cd LAGraph/build
//      cmake ..
//      make -j8
//
// Then run this demo with an input matrix.  For example:
//
//      ./experimental/benchmark/multisourcebfs_demo < ../data/west0067.mtx
//      ./experimental/benchmark/multisourcebfs_demo < ../data/karate.mtx
//      set up for raid graphs:
//      ./experimental/benchmark/multisourcebfs_demo < /raid/matrices/com-Youtube/com-Youtube.mtx
//
#include "LAGraphX.h"
#include "LG_internal.h"

#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&Y) ;                             \
    GrB_free (&A) ;                             \
    GrB_free (&parent) ;                        \
    GrB_free (&level) ;                         \
    GrB_free (&check_parent) ;                  \
    GrB_free (&check_level) ;                   \
    GrB_free (&SourceNodes) ;                   \
    LAGraph_Delete (&G, msg) ;                  \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Matrix Y = NULL, A = NULL ;
    GrB_Matrix level = NULL ;
    GrB_Matrix parent = NULL ;
    GrB_Matrix check_level = NULL ;
    GrB_Matrix check_parent = NULL ;
    GrB_Vector SourceNodes = NULL ;
    int nSources = 34;

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
    // try the LAGraph_MultiSourceBFS algorithm
    //--------------------------------------------------------------------------

    printf ("\n==========================Set up for BFS\n") ;
    // set up
    GRB_TRY (GrB_Vector_new (&SourceNodes, GrB_INT32, nSources)) ;
    printf ("\n==========================Intermediate Print\n") ;
    for (int i = 0; i < nSources; i++){
        GRB_TRY (GrB_Vector_setElement (SourceNodes, i, i)) ; 
    }
    

    printf ("\n==========================Running BFS\n") ;
    t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, &parent, G, SourceNodes, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_MultiSourceBFS: %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // check the results (make sure outputs match)
    //--------------------------------------------------------------------------

    bool levelisequal, parentisequal = false;
    GrB_Index n;
    GRB_TRY (GrB_Matrix_ncols (&n, level)) ;
    GrB_Type int_type = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    GRB_TRY (GrB_Matrix_new (&check_level, int_type, nSources, n)) ;
    GRB_TRY (GrB_Matrix_new (&check_parent, int_type, nSources, n)) ;
    t = LAGraph_WallClockTime ( ) ;
    for (int i = 0; i < nSources; i++) {
        GrB_Vector tempLevel, tempParent = NULL;
        GrB_Index src;
        GRB_TRY (GrB_Vector_extractElement (&src, SourceNodes, i)) ;
        LAGRAPH_TRY (LAGr_BreadthFirstSearch (&tempLevel,&tempParent, G, src, msg)) ;
        
        LAGRAPH_TRY (GrB_assign(check_level, NULL, NULL, tempLevel, i, GrB_ALL, n, GrB_NULL)) ;
        LAGRAPH_TRY (GrB_assign(check_parent, NULL, NULL, tempParent, i, GrB_ALL, n, GrB_NULL)) ;
    }
    
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to run equivalent regular BFSs:       %g sec\n", t) ;
    LAGRAPH_TRY (LAGraph_Matrix_IsEqual (&levelisequal, check_level, level, msg)) ;
    LAGRAPH_TRY (LAGraph_Matrix_IsEqual (&parentisequal, check_parent, parent, msg)) ;
    if (levelisequal ) // parents can be different in BFS
    {
        printf ("Test passed.\n") ;
    }
    else
    {
        printf ("Test failure!\n") ;
    }

    //--------------------------------------------------------------------------
    // print the results 
    //--------------------------------------------------------------------------

    printf ("\n===============================The result matrix level:\n") ;
    LAGRAPH_TRY (LAGraph_Matrix_Print (level, LAGraph_SHORT, stdout, msg)) ;
    printf ("\n===============================The result matrix parent:\n") ;
    LAGRAPH_TRY (LAGraph_Matrix_Print (parent, LAGraph_SHORT, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // free everything and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}

