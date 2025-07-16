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
#include "LG_Xtest.h"

void LG_rcc_iseq_demo(bool *z, const double *x, const double *y)
{
    (*z) = (isnan(*x) && isnan(*y)) ||*x == *y ;
}
#define ISEQ \
"   void LG_rcc_iseq_demo(bool *z, const double *x, const double *y) \n"\
"   {                                                           \n"\
"       (*z) = (isnan(*x) && isnan(*y)) || *x == *y ;           \n"\
"   }"
// LG_FREE_ALL is required by LG_TRY
#undef  LG_FREE_ALL
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&rcc1) ;                          \
    GrB_free (&rcc2) ;                          \
    GrB_free (&iseqFP) ;                        \
    LAGraph_Delete (&G, msg) ;                  \
}
#define SINGLERCC 1
int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;        // for error messages from LAGraph
    LAGraph_Graph G = NULL ;
    GrB_Vector rcc1 = NULL, rcc2 = NULL ;
    GrB_BinaryOp iseqFP = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = true ;               // set true for diagnostic outputs
    demo_init (burble) ;
    GRB_TRY (GxB_BinaryOp_new (
        &iseqFP, (GxB_binary_function) LG_rcc_iseq_demo, 
        GrB_BOOL, GrB_FP64, GrB_FP64, "LG_rcc_iseq_demo", ISEQ)) ;
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
    // try the LAGraph_RCC algorithm
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Cached_OutDegree (G, msg)) ;
    printf ("\n========================== Start RCC ==========================\n") ;
    t = LAGraph_WallClockTime ( ) ;
    LG_TRY (LAGraph_RichClubCoefficient (&rcc1, G, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_RichClubCoefficient: %g sec\n", t) ;
    
    #if SINGLERCC
    int result;
    LG_TRY (LAGraph_Cached_OutDegree (G, msg)) ;
    printf ("\n========================== Start RCC ==========================\n") ;
    t = LAGraph_WallClockTime ( ) ;
    result = LG_check_rcc (&rcc2, G, msg) ;
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time for LAGraph_RichClubCoefficient: %g sec\n", t) ;
    #endif

    //--------------------------------------------------------------------------
    // check the results (make sure Y is a copy of G->A)
    //--------------------------------------------------------------------------

    t = LAGraph_WallClockTime ( ) ;
    #if SINGLERCC
    if(result == GrB_SUCCESS)
    {
        bool flag;
        LG_TRY (LAGraph_Vector_IsEqualOp(&flag, rcc1, rcc2, iseqFP, msg)) ;
        if (flag)
            printf("TEST PASSED\n") ;
        else
            printf("TEST FAILED\n") ;
    }
    else
        printf("Test indeterminate. Single Thread exited with %d\n", result) ;
    #endif
    t = LAGraph_WallClockTime ( ) - t ;
    printf ("Time to check results:       %g sec\n", t) ;

    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf ("\n===============================The result matrix Y:\n") ;
    GRB_TRY (GxB_Vector_fprint (rcc1, "rcc", GxB_SHORT, stdout));

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
