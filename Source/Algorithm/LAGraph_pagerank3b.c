//------------------------------------------------------------------------------
// LAGraph_pagerank3b: pagerank using a real semiring
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

// LAGraph_pagerank3b: Alternative PageRank implementation using a real
// semiring.
//
// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL {       \
    GrB_free(&transpose_desc);   \
    GrB_free(&invmask_desc);     \
    GrB_free(&A);                \
    GrB_free(&G);                \
    GrB_free(&grb_d_out);            \
    GrB_free(&importance_vec);   \
    GrB_free(&grb_pr);               \
};

// uncomment this to see the intermidiate resluts; lots of prints!!
//#undef NDEBUG

// uncomment this to see the  timing info
#define PRINT_TIMING_INFO

    GrB_Info LAGraph_pagerank3b // PageRank definition
(
 GrB_Vector *result,    // output: array of LAGraph_PageRank structs
 GrB_Matrix A,          // binary input graph, not modified
 float damping_factor, // damping factor
 unsigned long itermax, // maximum number of iterations
 int* iters             // output: number of iterations taken
 )
{
    GrB_Info info;
    GrB_Index n;

    GrB_Descriptor invmask_desc;
    GrB_Descriptor transpose_desc;
    GrB_Vector grb_d_out;

#ifdef PRINT_TIMING_INFO
    // start the timer
    double tic [2] ;
    LAGraph_tic (tic) ;
#endif


    GrB_Vector importance_vec = NULL ;
    GrB_Vector grb_pr = NULL;

    GrB_Matrix G; // a dense row of zeros zeroes(1,nc)
    GrB_Index nc; //number of columnns

    LAGRAPH_OK(GrB_Matrix_ncols(&nc, A));

    LAGRAPH_OK(GrB_Matrix_nrows(&n, A));
    GrB_Index nvals;
    LAGRAPH_OK(GrB_Matrix_nvals(&nvals, A));


    LAGRAPH_OK(GrB_Matrix_new (&G, GrB_FP32, n, nc));

    // G is zeros in last row
    for (GrB_Index c = 0; c < nc; c++){
        LAGRAPH_OK(GrB_Matrix_setElement (G, 0.0, nc-1, c));
    }
#ifndef NDEBUG
    int print_size = 5;   //number of entries get printed
    print_size = (print_size > n)? n : print_size;
    //  GxB_print (G, 3) ;
#endif

    // A += G;
    LAGRAPH_OK(GrB_eWiseAdd (A, NULL, NULL, GrB_PLUS_FP32, A, G, NULL));


    LAGRAPH_OK(GxB_set (A, GxB_FORMAT, GxB_BY_COL));

#ifndef NDEBUG
    // GxB_print (A, 3) ;
#endif


    // Create complement descriptor
    LAGRAPH_OK(GrB_Descriptor_new(&invmask_desc));
    LAGRAPH_OK(GrB_Descriptor_set(invmask_desc, GrB_MASK, GrB_SCMP));

    // Create transpose descriptor
    LAGRAPH_OK(GrB_Descriptor_new(&transpose_desc));
    LAGRAPH_OK(GrB_Descriptor_set(transpose_desc, GrB_INP0, GrB_TRAN));
    LAGRAPH_OK(GrB_Descriptor_set(transpose_desc, GrB_OUTP, GrB_REPLACE));

    // Matrix A row sum
    // Stores the outbound degrees of all vertices
    LAGRAPH_OK(GrB_Vector_new(&grb_d_out, GrB_UINT64, n));
    LAGRAPH_OK(GrB_reduce( grb_d_out, NULL, NULL, GxB_PLUS_UINT64_MONOID,
                A, NULL ));

#ifndef NDEBUG
    GxB_print (grb_d_out, 1) ;
    // GxB_print (A, 3) ;
#endif

    // Iteration
    // Initialize PR vector
    LAGRAPH_OK(GrB_Vector_new(&grb_pr, GrB_FP32, n));

    LAGRAPH_OK(GrB_Vector_new(&importance_vec, GrB_FP32, n));

    // Teleport value
    const float teleport = (1 - damping_factor) / n;

    float  tol = 1e-4;
    float  rdiff = 1 ;       // first iteration is always done

    GrB_Type type = GrB_FP32 ;
    GrB_Index *dI = NULL ;
    unsigned long int *d_sp= NULL ;
    GrB_Index d_nvals;
    GrB_Index d_n;

    // d_sp <-----   grb_d_out || export
    LAGRAPH_OK (GxB_Vector_export (&grb_d_out, &type, &d_n, &d_nvals, &dI,
                (void **) (&d_sp),   NULL)) ;

    // dens d_out
    long int *d_out = (long int*) calloc(n, sizeof(long int));

    int nthreads = LAGraph_get_nthreads ( ) ;
    nthreads = LAGRAPH_MIN (n , nthreads) ;
    nthreads = LAGRAPH_MAX (nthreads, 1) ;

    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i = 0 ; i < d_nvals; i++){
        GrB_Index ind = (GrB_Index) dI[i];
        d_out [ind] = d_sp [i];
    }

    free (d_sp);
    free (dI);

#ifndef NDEBUG
    for (int i = 0 ; i < print_size; i++){
        printf("d_out [%d]=%ld\n", i, d_out [i]);
    }
#endif



    // initializing pr
    float *pr = (float *) malloc (n*sizeof(float));     
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i = 0; i < n ; i++){
        pr [i] = 1.0/n;
   }
#ifndef NDEBUG
    for (int i = 0 ; i < print_size ; i++){
        printf("pr[%d]=%f\n", i, pr [i]);
    }
#endif
 
    float *oldpr = (float *) malloc (n*sizeof(float));     

    //initailze the dense indices
    GrB_Index *I = LAGraph_malloc(n, sizeof(GrB_Index));
    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (GrB_Index j = 0; j < n; j++){
        I[j] = j;
    }

#ifdef PRINT_TIMING_INFO
    // stop the timer
    double t1 = LAGraph_toc (tic); 
    printf ("\ninitialization time: %12.6e (sec)\n",t1);
    LAGraph_tic (tic);
#endif



    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++) {
        // oldpr = pr; deep copy
        //GrB_Vector_dup(&oldpr, pr);
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0; i < n ; i++){
            oldpr [i] = pr [i];
        }

        // Importance calculation
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0 ; i < n; i++){
            if (d_out [i] != 0){
                pr [i] = damping_factor * pr [i] / d_out [i];
            }
            else{ 
                pr [i] = 0; 
            }
        }


#ifndef NDEBUG
        for (int i = 0 ; i < print_size; i++){
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
        // importance_vec *= importance_vec * A'?
        LAGRAPH_OK(GrB_mxv( importance_vec, NULL, NULL, GxB_PLUS_TIMES_FP32,
                    A, importance_vec, transpose_desc ));
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
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (int i = 0 ; i < n; i++){
            pr [i] +=  teleport; 

        }

#ifndef NDEBUG
        for (int i = 0 ; i < print_size; i++){
            printf (" pr [%d] = %f\n",  i,  pr [i]);
        }
#endif

        //----------------------------------------------------------------------
        // rdiff = sum ((oldpr-pr).^2)
        //----------------------------------------------------------------------
        rdiff = 0;
        // norm (oldpr pr, 1)
        #pragma omp parallel for num_threads(nthreads) reduction(+:rdiff)
        for (int i = 0 ; i < n; i++){
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
    printf ("compuatatin time: %12.6e (sec) ratio (comp/init): %f\n\n",
            t2, t2/t1);
#endif



    GrB_Index *prI = LAGraph_malloc(n, sizeof(GrB_Index));

    // grb_pr<-----  pr || import back
    LAGRAPH_OK (GxB_Vector_import (&grb_pr, GrB_FP32, n, n, &I, 
                (void **) (&pr),   NULL)) ;
    (*result) = grb_pr;

    free(I);
    free (oldpr);
    return (GrB_SUCCESS);
}
