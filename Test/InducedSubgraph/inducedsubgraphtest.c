//------------------------------------------------------------------------------
// LAGraph/Test/LCC/extractsubmatrixtest.c: test program for LAGraph_extractsubmatrix
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

// Contributed by Tim Davis, Texas A&M and Gabor Szarnyas, BME

// Usage: inducedsubgraphtest can be used with binary input graphs
//
// inducedsubgraphtest binarymatrixfile.grb

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL                            \
{                                                   \
    GrB_free (&A) ;                                 \
    GrB_free (&C) ;                                 \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Matrix A = NULL, C = NULL ;

    LAGraph_init ( ) ;
    int nthreads_max = LAGraph_get_nthreads ( ) ;
    if (nthreads_max == 0) nthreads_max = 1 ;

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    if (argc < 1)
    {
        printf ("Usage: inducedsubgraphtest binarymatrixfile.grb\n") ;
        return (GrB_INVALID_VALUE) ;
    }

    LAGRAPH_OK (LAGraph_binread(&A, argv[1])) ;


    //--------------------------------------------------------------------------
    // extract induced subgraph
    // - method 1: multiply matrix from left and right with diagm(nodes)
    // - method 2: use select operator
    //--------------------------------------------------------------------------

    #define NTRIALS 7
    int nthread_list [NTRIALS] = { 1, 2, 4, 8, 16, 32, 64 } ;

    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        int nthreads = nthread_list [trial] ;
        if (nthreads > nthreads_max) break ;
        LAGraph_set_nthreads (nthreads) ;

        // ignore the sanitize time;  assume the user could have provided an
        // input graph that is already binary with no self-edges

        GrB_Matrix C;
        GrB_Index n;
        GrB_Matrix_nrows(&n, A);

        // select every other vertex in the graph
        GrB_Index nv = n / 2;
        GrB_Index* V = LAGraph_malloc(nv, sizeof(GrB_Index));
        for (int k = 0; k < nv; k++) {
            V[k] = 2*k;
        }

        double tic [2] ;
        LAGraph_tic (tic) ;
        LAGRAPH_OK (LAGraph_inducedsubgraph(&C, A, V, nv, true)) ;
//        LAGRAPH_OK (LAGraph_inducedsubgraph(&C, A, V, nv, false)) ;
        double time = LAGraph_toc (tic) ;
        printf("Time elapsed: %10.2f seconds, %d threads\n", time, nthreads) ;

        //LAGRAPH_OK ( GxB_print(C, GxB_SHORT)) ;

        LAGRAPH_FREE(C);
    }

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}
