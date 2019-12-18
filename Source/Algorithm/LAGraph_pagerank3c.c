//------------------------------------------------------------------------------
// LAGraph_pagerank3c: pagerank using a real semiring
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

// LAGraph_pagerank3c: Alternative PageRank implementation using a real
// semiring.
//
// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619

// This algorithm assumes the graph has no nodes with no out-going edges.  In
// terms of the adjacency matrix, it assumes there are no rows in A that have
// no entries.

// For fastest results, the input matrix should stored in GxB_BY_COL format.

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL        \
{                               \
    LAGRAPH_FREE (I);           \
    LAGRAPH_FREE (d_out) ;      \
    GrB_free(&grb_d_out);       \
    GrB_free(&importance_vec);  \
    GrB_free(&grb_pr);          \
};

#undef NDEBUG
// comment this out to see the intermediate resluts; lots of prints!!
#define NDEBUG

// uncomment this to see the  timing info
#define PRINT_TIMING_INFO

GrB_Info LAGraph_pagerank3c // PageRank definition
(
 GrB_Vector *result,    // output: array of LAGraph_PageRank structs
 GrB_Matrix A,          // binary input graph, not modified
 float damping_factor,  // damping factor
 unsigned long itermax, // maximum number of iterations
 int* iters             // output: number of iterations taken
 )
{

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Info info;
    GrB_Index n;

    GrB_Vector grb_d_out = NULL ;
    GrB_Vector importance_vec = NULL ;
    GrB_Vector grb_pr = NULL;
    GrB_Index *I = NULL ;
    float *d_out = NULL ;

#ifdef PRINT_TIMING_INFO
    // start the timer
    double tic [2] ;
    LAGraph_tic (tic) ;
#endif

    GrB_Index ncols ; //number of columnns

    LAGRAPH_OK(GrB_Matrix_ncols(&ncols , A));
    LAGRAPH_OK(GrB_Matrix_nrows(&n, A));
    GrB_Index nvals;
    LAGRAPH_OK(GrB_Matrix_nvals(&nvals, A));

    if (ncols  != n)
    {
        return (GrB_DIMENSION_MISMATCH) ;
    }

    // Matrix A row sum
    // Stores the outbound degrees of all vertices
    LAGRAPH_OK(GrB_Vector_new(&grb_d_out, GrB_FP32, n));
    LAGRAPH_OK(GrB_reduce( grb_d_out, NULL, NULL, GxB_PLUS_FP32_MONOID,
                A, NULL ));

#ifndef NDEBUG
    GxB_print (grb_d_out, 1) ;
    // GxB_print (A, 3) ;
#endif

    // Teleport value
    const float teleport = (1 - damping_factor) / n;

    float  tol = 1e-4;
    float  rdiff = 1 ;       // first iteration is always done

    GrB_Type type = GrB_FP32 ;
    GrB_Index d_nvals;
    GrB_Index d_n;

    // d_out <-----   grb_d_out || export
    LAGRAPH_OK (GxB_Vector_export (&grb_d_out, &type, &d_n, &d_nvals, &I,
                (void **) (&d_out),   NULL)) ;

    if (d_nvals < n)
    {
        LAGRAPH_ERROR ("Matrix has dangling nodes!", GrB_INVALID_VALUE) ;
    }

    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (n , nthreads) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

#ifndef NDEBUG
    for (int i = 0 ; i < n; i++){
        printf("d_out [%d]=%ld\n", i, d_out [i]);
    }
#endif

    //--------------------------------------------------------------------------
    // compute the pagerank
    //--------------------------------------------------------------------------

#ifdef PRINT_TIMING_INFO
    // stop the timer
    double t1 = LAGraph_toc (tic); 
    printf ("\ninitialization time: %12.6e (sec)\n",t1);
    LAGraph_tic (tic);
#endif

    // initializing pr
    float *pr = (float *) malloc (n*sizeof(float));     
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i = 0; i < n ; i++){
        pr [i] = 1.0/n;
   }
#ifndef NDEBUG
    for (int i = 0 ; i < n ; i++){
        printf("pr[%d]=%f\n", i, pr [i]);
    }
#endif
 
    float *oldpr = (float *) malloc (n*sizeof(float));     

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {

        // Importance calculation
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0 ; i < n; i++){
            oldpr [i] = pr [i];
            // assumes d_out [i] > 0
            pr [i] = damping_factor * pr [i] / d_out [i];
        }

#ifndef NDEBUG
        for (int i = 0 ; i < n; i++){
            printf (" pr [%d] = %f\n",  i,  pr [i]);
        }
#endif
        // importance_vec <-----  pr
        LAGRAPH_OK (GxB_Vector_import (&importance_vec, GrB_FP32, n, n, &I, 
                    (void **) (&pr),   NULL)) ;

#ifndef NDEBUG
        printf ("after importance_vec import\n");
        GxB_print (importance_vec, 2) ;
#endif

        // Calculate total PR of all inbound vertices
        // importance_vec = A' * importance_vec
        LAGRAPH_OK(GrB_mxv( importance_vec, NULL, NULL,
            GxB_PLUS_SECOND_FP32, A, importance_vec, LAGraph_desc_toor ));

#ifndef NDEBUG
        printf ("==============2\n");
        printf ("after mxv\n");
        GxB_print (importance_vec, 1) ;
#endif

        GrB_Index nvals_exp;
        // pr <-----  importance_vec 
        GrB_Type ivtype;
        LAGRAPH_OK (GxB_Vector_export (&importance_vec, &ivtype, &n, &nvals_exp,
                    &I,  (void **) (&pr),   NULL)) ;
        // assert (nvals_exp == n );


        // PageRank summarization
        // Add teleport, importance_vec, and dangling_vec components together
        // pr = (1-df)/n + pr
        // rdiff = sum (abs(oldpr-pr))

        rdiff = 0;
        // norm (oldpr pr, 1)
        #pragma omp parallel for num_threads(nthreads) schedule(static) \
            reduction(+:rdiff)
        for (int i = 0 ; i < n; i++)
        {
            pr [i] += teleport; 
            float d = (oldpr [i] - pr [i]); 
            d = (d > 0 ? d : -d); //abs(d)
            rdiff += d;
        }

#ifndef NDEBUG
        printf("---------------------------iters %d  rdiff=%f\n",*iters, rdiff);
#endif
    }

#ifdef PRINT_TIMING_INFO
    // stop the timer
    double t2 = LAGraph_toc (tic); 
    printf ("computation time: %12.6e (sec) ratio (comp/init): %f\n\n",
            t2, t2/t1);
#endif

    // grb_pr<-----  pr || import back
    LAGRAPH_OK (GxB_Vector_import (&grb_pr, GrB_FP32, n, n, &I, 
                (void **) (&pr),   NULL)) ;
    (*result) = grb_pr;
    grb_pr = NULL ;

    LAGRAPH_FREE (oldpr);
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS);
}

