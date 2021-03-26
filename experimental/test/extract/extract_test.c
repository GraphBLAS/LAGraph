//------------------------------------------------------------------------------
// LAGraph/Test/extract/extract_test.c: test program for GrB*extractElement
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

// Contributed by Tim Davis, Texas A&M

// Usage:  extract_test

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include "LAGraph.h"

#define LAGRAPH_FREE_ALL    \
    GrB_free (&X) ;

#define CHECK(ok)                                                   \
    if (!(ok))                                                      \
    {                                                               \
        fprintf (stderr, "fail: %s %d\n", __FILE__, __LINE__) ;     \
        printf ("fail: %s %d\n", __FILE__, __LINE__) ;              \
        abort ( ) ;                                                 \
    }

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Vector X = NULL ;
    LAGraph_init ( ) ;         // start LAGraph

    //--------------------------------------------------------------------------
    // construct a nearly-dense vector
    //--------------------------------------------------------------------------

    double tic [2] ;
    LAGraph_tic (tic) ;

    GrB_Index n = 64 * 1024 * 1024 ;
    printf ("extract test: n = %lu\n", n) ;
    LAGr_Vector_new (&X, GrB_UINT64, n) ;
    for (uint64_t k = 1 ; k < n ; k++)
    {
        LAGr_Vector_setElement (X, k, k) ;
    }

    double t = LAGraph_toc (tic) ;
    printf ("set     time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    LAGraph_tic (tic) ;
    GrB_Index ignore ;
    LAGr_Vector_nvals (&ignore, X) ;
    t = LAGraph_toc (tic) ;
    printf ("wait    time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    LAGraph_tic (tic) ;
    GxB_print (X, GxB_SHORT) ;
    t = LAGraph_toc (tic) ;
    printf ("check   time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    //--------------------------------------------------------------------------
    // test binary searches
    //--------------------------------------------------------------------------

    uint64_t x ;

    for (uint64_t k = 1 ; k < n ; k++)
    {
        LAGr_Vector_extractElement (&x, X, k) ;
        CHECK (x == k) ;
    }

    t = LAGraph_toc (tic) ;
    printf ("extract time %12.6f n/sec %12.6f million\n", t, 1e-6 * ((double) n-1) / t) ;

    info = GrB_Vector_extractElement (&x, X, 0) ;
    CHECK (info == GrB_NO_VALUE) ;

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
    return (GrB_SUCCESS) ;
}

