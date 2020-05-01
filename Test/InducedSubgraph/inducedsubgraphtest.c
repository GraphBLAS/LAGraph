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

    GrB_Index n;
    GrB_Matrix_nrows(&n, A);

    for (GrB_Index k = 2; k <= 2; k++) {
        // select every kth vertex in the graph
        GrB_Index nv = n / k;

        //--------------------------------------------------------------------------
        // extract induced subgraph
        // - method 1: multiply matrix from left and right with diagm(nodes)
        // - method 2: use select operator
        //--------------------------------------------------------------------------

        #define NTRIALS 7
        int nthread_list[NTRIALS] = {1, 2, 4, 8, 16, 32, 64};

        for (int trial = 0; trial < NTRIALS; trial++) {
            int nthreads = nthread_list[trial];
            if (nthreads > nthreads_max) break;
            LAGraph_set_nthreads(nthreads);

            GrB_Matrix C = NULL;

            double tic[2];
            LAGraph_tic(tic);
            GrB_Index *Vsparse = LAGraph_malloc(nv, sizeof(GrB_Index));
            for (int i = 0; i < nv; i++) {
                Vsparse[i] = k * i;
            }
            LAGRAPH_OK(LAGraph_Matrix_extract_keep_dimensions(&C, A, Vsparse, NULL, nv));
            LAGRAPH_FREE(Vsparse);
            double time = LAGraph_toc(tic);
            printf("Vsparse\t%ld\t%d\t%.2f\n", k, nthreads, time);

            LAGRAPH_FREE(C);
        }

        for (int trial = 0; trial < NTRIALS; trial++) {
            int nthreads = nthread_list[trial];
            if (nthreads > nthreads_max) break;
            LAGraph_set_nthreads(nthreads);

            GrB_Matrix C = NULL;

            double tic[2];
            LAGraph_tic(tic);
            bool *Vdense = LAGraph_calloc(nv, sizeof(bool));
            for (int i = 0; i < nv; i++) {
                Vdense[i] = true;
            }
            LAGRAPH_OK(LAGraph_Matrix_extract_keep_dimensions(&C, A, NULL, Vdense, nv));
            LAGRAPH_FREE(Vdense);
            double time = LAGraph_toc(tic);
            printf("Vdense\t%ld\t%d\t%.2f\n", k, nthreads, time);

            LAGRAPH_FREE(C);
        }
    }

    LAGRAPH_FREE_ALL ;
    LAGraph_finalize ( ) ;
}
