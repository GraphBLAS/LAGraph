//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/hdip_fielder_demo.c: a simple demo
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_internal.h"

float difference (GrB_Vector centrality, double *matlab_result) ;

float difference (GrB_Vector centrality, double *matlab_result)
{
    GrB_Vector diff = NULL, cmatlab = NULL ;
    GrB_Index n = 0 ;
    GrB_Vector_size (&n, centrality) ;
    GrB_Vector_new (&cmatlab, GrB_FP32, n) ;
    for (int i = 0 ; i < n ; i++)
    {
        GrB_Vector_setElement_FP64 (cmatlab, matlab_result [i], i) ;
    }
    // diff = max (abs (cmatlab - centrality))
    GrB_Vector_new (&diff, GrB_FP32, n) ;
    GrB_eWiseAdd (diff, NULL, NULL, GrB_MINUS_FP32, cmatlab, centrality,
        NULL) ;
    GrB_apply (diff, NULL, NULL, GrB_ABS_FP32, diff, NULL) ;
    float err = 0 ;
    GrB_reduce (&err, NULL, GrB_MAX_MONOID_FP32, diff, NULL) ;
    GrB_free (&diff) ;
    GrB_free (&cmatlab) ;
    return (err) ;
}

// Bucky test data
double bucky_fiedler [60] = {
	-0.2236,
	-0.2071,
	-0.1804,
	-0.1804,
	-0.2071,
	-0.2022,
	-0.1669,
	-0.1098,
	-0.1098,
	-0.1669,
	-0.1669,
	-0.1481,
	-0.0744,
	-0.0477,
	-0.1049,
	-0.1098,
	-0.0744,
	 0.0094,
	 0.0259,
   	-0.0477,
   	-0.1098,
   	-0.0477,
   	 0.0259,
   	 0.0094,
   	-0.0744,
   	-0.1669,
   	-0.1049,
   	-0.0477,
   	-0.0744,
   	-0.1481,
   	 0.1481,
   	 0.0745,
   	 0.0477,
   	 0.1049,
   	 0.1669,
    	 0.0745,
   	-0.0094,
   	-0.0259,
   	 0.0477,
   	 0.1098,
   	 0.0477,
   	-0.0259,
   	-0.0094,
   	 0.0745,
   	 0.1098,
   	 0.1049,
   	 0.0477,
   	 0.0745,
   	 0.1481,
   	 0.1669,
   	 0.1669,
   	 0.1098,
   	 0.1098,
   	 0.1669,
   	 0.2022,
   	 0.2071,
   	 0.1804,
   	 0.1804,
   	 0.2071,
   	 0.2236} ;

// Karate test data
double karate_fiedler [34] = {
        -0.3561,
        -0.1036,
        -0.0156,
        -0.1243,
        -0.2280,
        -0.2097,
        -0.2097,
        -0.1224,
         0.0163,
         0.1108,
        -0.2280,
        -0.2463,
        -0.1853,
        -0.0725,
         0.1900,
         0.1900,
        -0.1548,
        -0.1749,
         0.1900,
        -0.0741,
         0.1900,
        -0.1749,
         0.1900,
         0.1792,
         0.1703,
         0.1794,
         0.2155,
         0.1428,
         0.1002,
         0.1937,
         0.0732,
         0.0790,
         0.1427,
         0.1274} ;

double west0067_fiedler [67] = {
   -0.7918,
   -0.0506,
   -0.0329,
   -0.0366,
   -0.1569,
   -0.1608,
   -0.1776,
   -0.1747,
   -0.1529,
   -0.0391,
   -0.0320,
   -0.0046,
   -0.1138,
   -0.0140,
   -0.0314,
   -0.0114,
   -0.0066,
   -0.0862,
    0.0245,
   -0.0117,
    0.0232,
    0.0338,
    0.0052,
    0.0185,
   -0.0731,
   -0.0520,
   -0.0602,
   -0.0711,
   -0.0623,
    0.0528,
   -0.0016,
    0.0447,
    0.0566,
    0.0444,
    0.0610,
    0.0220,
   -0.0008,
    0.0171,
    0.0305,
    0.0519,
    0.0414,
    0.0491,
    0.0482,
    0.0912,
    0.0660,
    0.1074,
    0.1016,
    0.1078,
    0.0683,
    0.0871,
    0.0777,
    0.0839,
    0.0901,
    0.1092,
    0.0850,
    0.0752,
   -0.0019,
    0.0239,
    0.0442,
    0.0772,
   -0.0179,
    0.0770,
    0.1072,
    0.0342,
    0.0762,
    0.1115,
    0.1000} ;

// LG_FREE_ALL is required by LG_TRY
#undef LG_FREE_ALL
#define LG_FREE_ALL              \
    {                            \
        GrB_free(&Y);            \
        GrB_free(&x);            \
        GrB_free(&u);            \
        GrB_free(&steper);       \
        GrB_free(&indiag);       \
        LAGraph_Delete(&G, msg); \
        LAGraph_Finalize(msg);   \
    }

int main(int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // startup LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg[LAGRAPH_MSG_LEN]; // for error messages from LAGraph
    LAGraph_Graph G = NULL;
    GrB_Matrix Y = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix indiag = NULL;
    GrB_Vector x = NULL;
    // Additional variables and modifications needed to test MYPCG2
    GrB_Vector steper = NULL;
    GrB_Vector u = NULL; // a vector of size nrowsLap, filled with 1.
    // set u[0] = 1+sqrt(nrowsLap)
    // Additional variables needed to test Hdip
    GrB_Vector iters_handle = NULL;
    float lambda_result = 0;
    GrB_Vector x_handle = NULL;
    GrB_Vector kmax = NULL;

    // start GraphBLAS and LAGraph
    bool burble = false; // set true for diagnostic outputs
    demo_init(burble);

    //--------------------------------------------------------------------------
    // read in the graph: this method is defined in LAGraph_demo.h
    //--------------------------------------------------------------------------

    // readproblem can read in a file in Matrix Market format, or in a binary
    // format created by binwrite (see LAGraph_demo.h, or the main program,
    // mtx2bin_demo).

    double t = LAGraph_WallClockTime();
    char *matrix_name = (argc > 1) ? argv[1] : "stdin";
    
    LG_TRY(readproblem(
        &G,           // the graph that is read from stdin or a file
        NULL,         // source nodes (none, if NULL)
        true,         // make the graph undirected, if true
        true,         // remove self-edges, if true
        false,        // return G->A as structural, if true,
        NULL,         // prefered GrB_Type of G->A; null if no preference
        false,        // ensure all entries are positive, if true
        argc, argv)); // input to this main program
    t = LAGraph_WallClockTime() - t;
    printf("Time to read the graph:      %g sec\n", t);

    printf("\n==========================The input graph matrix G:\n");
    LG_TRY(LAGraph_Graph_Print(G, LAGraph_SHORT, stdout, msg));

    //--------------------------------------------------------------------------
    // try the LAGraph_HelloWorld "algorithm"
    //--------------------------------------------------------------------------
    GrB_Index n, nvals;
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));
    GRB_TRY(GrB_Matrix_nvals(&nvals, G->A));
    GRB_TRY(GrB_Matrix_new(&A, GrB_FP32, n, n));
    GRB_TRY(GrB_assign(A, G->A, NULL, (double)1,
                       GrB_ALL, n, GrB_ALL, n, GrB_DESC_S));
    GrB_free(&(G->A));
    G->A = A;

    //--------------------------------------------------------------------------
    // try the LAGraph_Laplacian "algorithm"
    //--------------------------------------------------------------------------

    // Variables needed to test Laplacian
    float infnorm;
    LG_TRY(LAGraph_Laplacian(&Y, &infnorm, A, msg));

    printf("\n===========================The laplacian matrix: \n");
    LAGRAPH_TRY(LAGraph_Matrix_Print(Y, LAGraph_SHORT, stdout, msg));

    //--------------------------------------------------------------------------
    // try the LAGraph_mypcg2 "algorithm"
    //--------------------------------------------------------------------------

    GrB_Index k;
    GrB_Index nrows;
    float nrowsLap; // number of rows of laplacian matrix
    float alpha;
    GRB_TRY(GrB_Matrix_nrows(&nrows, Y));

    nrowsLap = (float)n;

    GRB_TRY(GrB_Vector_new(&u, GrB_FP32, n));
    // u = all ones vector
    GRB_TRY(GrB_assign(u, NULL, NULL, 1, GrB_ALL, n, NULL));
    // u [0] = 1+sqrt(n)
    GRB_TRY(GrB_Vector_setElement_FP32(u, 1 + sqrt(nrowsLap), 0));

    alpha = nrowsLap + sqrt(nrowsLap);

    GRB_TRY(GrB_Matrix_new(&indiag, GrB_FP32, n, n));
    GRB_TRY(GrB_select(indiag, NULL, NULL, GrB_DIAG, Y, 0, NULL));
    GRB_TRY(GrB_apply(indiag, NULL, NULL, GrB_MINV_FP32, indiag, NULL));
    
    GRB_TRY(GrB_Vector_new(&x, GrB_FP32, n));
    GRB_TRY(GrB_assign(x, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Vector_setElement_FP32(x, 0, 0));

    t = LAGraph_WallClockTime( );
    int result = LAGraph_mypcg2(&steper, &k, Y, u, alpha, indiag, x, .000001, 50, msg);
    t = LAGraph_WallClockTime( ) - t;
    printf("Time for LAGraph_mypcg2: %g sec\n", t);

    //--------------------------------------------------------------------------
    // try the LAGraph_Hdip_Fiedler "algorithm"
    //--------------------------------------------------------------------------
    //Set kmax = [20,50]
    GRB_TRY(GrB_Vector_new(&kmax, GrB_FP32, 2));
    GRB_TRY(GrB_Vector_setElement_FP32(kmax, 20, 0));
    GRB_TRY(GrB_Vector_setElement_FP32(kmax, 50, 1));

    t = LAGraph_WallClockTime( );
    LAGraph_Hdip_Fiedler(&iters_handle, &lambda_result, &x_handle, Y, infnorm, kmax, 0.000001, 0.000001, msg);
    t = LAGraph_WallClockTime( ) - t;
    printf("Time for LAGraph_Hdip_Fiedler: %g sec\n", t);

    //--------------------------------------------------------------------------
    // check the results (Make sure imported graph is correct for karate or bucky)
    //--------------------------------------------------------------------------
    
    //Testing karate
    t = LAGraph_WallClockTime();
    float err = difference (x_handle, karate_fiedler) ;
    printf("\n=============================Testing karate x vector:\n");
    t = LAGraph_WallClockTime() - t;
    printf("Time to check results:       %g sec\n", t);
    if (err < 1e-4)
    {
        printf("Test passed.\n");
    }
    else
    {
        printf("Test failure!\n");
    }

    printf("\n=============================Testing karate lambda:\n");
    err = fabs(lambda_result - 1.3297);
    if (err < 1e-4)
    {
        printf("Test passed.\n");
    }
    else
    {
        printf("Test failure!\n");
    }
    
    //Testing Bucky
    t = LAGraph_WallClockTime();
    err = difference (x_handle, bucky_fiedler) ;
    printf("\n=============================Testing bucky x vector:\n");
    t = LAGraph_WallClockTime() - t;
    printf("Time to check results:       %g sec\n", t);
    if (err < 1e-4)
    {
        printf("Test passed.\n");
    }
    else
    {
        printf("Test failure!\n");
    }

    printf("\n=============================Testing bucky lambda:\n");
    err = fabs(lambda_result - 0.2434);
    if (err < 1e-4)
    {
        printf("Test passed.\n");
    }
    else
    {
        printf("Test failure!\n");
    }

    //Testing West0067
    t = LAGraph_WallClockTime();
    err = difference (x_handle, west0067_fiedler) ;
    printf("\n=============================Testing west0067 x vector:\n");
    t = LAGraph_WallClockTime() - t;
    printf("Time to check results:       %g sec\n", t);
    if (err < 1e-4)
    {
        printf("Test passed.\n");
    }
    else
    {
        printf("Test failure!\n");
    }

    printf("\n=============================Testing west0067 lambda:\n");
    err = fabs(lambda_result - 6.5586);
    if (err < 1e-4)
    {
        printf("Test passed.\n");
    }
    else
    {
        printf("Test failure!\n");
    }


    //--------------------------------------------------------------------------
    // print the results (Y is just a copy of G->A)
    //--------------------------------------------------------------------------

    printf("\n===============================The result vector x:\n");
    LAGraph_Vector_Print (x_handle, 3, stdout, msg) ;
    printf("\n===============================The lambda: %f\n", lambda_result);
    printf("\n===============================The iters: \n");
    LAGraph_Vector_Print (iters_handle, 3, stdout, msg) ;

    //--------------------------------------------------------------------------
    // free everyting and finish
    //--------------------------------------------------------------------------

    LG_FREE_ALL;
    printf("finalize\n");
    LG_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS);
}
