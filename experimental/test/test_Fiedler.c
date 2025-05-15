//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_Fiedler.c: test the HDIP Fiedler method
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A Davis, Texas A&M University

//------------------------------------------------------------------------------

#include <acutest.h>
#include <LAGraph_test.h>
#include "LG_Xtest.h"
#include "LG_internal.h"

typedef struct
{
    const char *name ;
    double lambda ;
    double *fiedler ;
}
matrix_info ;

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

const matrix_info files [ ] =
{
    { "bucky.mtx",    0.2434, bucky_fiedler }, 
    { "karate.mtx",   1.3297, karate_fiedler },
    { "west0067.mtx", 6.5586, west0067_fiedler },
    { "", 0, NULL },
} ;

void test_fiedler (void)
{
#if LAGRAPH_SUITESPARSE

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
    GrB_Vector iters = NULL;
    float lambda_result = 0;
    GrB_Vector fiedler_vector = NULL;
    GrB_Vector kmax = NULL;
    #define LEN 512
    char filename [LEN+1] ;
    GrB_Index n ;

    // start GraphBLAS and LAGraph
    bool burble = false; // set true for diagnostic outputs
    LAGraph_Init (msg) ;

    //--------------------------------------------------------------------------
    // read in the graphs and test them
    //--------------------------------------------------------------------------

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        printf ("\n %s: ==================================\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;
        fclose (f) ;

        // set all entries to 1
        OK (GrB_Matrix_nrows (&n, A)) ;
        OK (GrB_assign(A, A, NULL, (double)1,
                           GrB_ALL, n, GrB_ALL, n, GrB_DESC_S));

        // ensure A is symmetric and remove self edges
        OK (GrB_eWiseAdd (A, NULL, NULL, GrB_ONEB_FP64, A, A, GrB_DESC_T1)) ;
        OK (GrB_select (A, NULL, NULL, GrB_OFFDIAG, A, 0, NULL)) ;

        //----------------------------------------------------------------------
        // try the LAGraph_Laplacian algorithm
        //----------------------------------------------------------------------

        // Variables needed to test Laplacian
        float infnorm;
        OK (LAGraph_Laplacian(&Y, &infnorm, A, msg));

        printf("\n===========================The laplacian matrix: \n");
        OK (LAGraph_Matrix_Print(Y, LAGraph_SHORT, stdout, msg)) ;

        //----------------------------------------------------------------------
        // try the LAGraph_mypcg2 algorithm
        //----------------------------------------------------------------------

        GrB_Index kk;
        GrB_Index nrows;
        float nrowsLap; // number of rows of laplacian matrix
        float alpha;
        OK (GrB_Matrix_nrows(&nrows, Y));

        nrowsLap = (float)n;

        OK (GrB_Vector_new(&u, GrB_FP32, n));
        // u = all ones vector
        OK (GrB_assign(u, NULL, NULL, 1, GrB_ALL, n, NULL));
        // u [0] = 1+sqrt(n)
        OK (GrB_Vector_setElement_FP32(u, 1 + sqrt(nrowsLap), 0));

        alpha = nrowsLap + sqrt(nrowsLap);

        OK (GrB_Matrix_new(&indiag, GrB_FP32, n, n));
        OK (GrB_select(indiag, NULL, NULL, GrB_DIAG, Y, 0, NULL));
        OK (GrB_apply(indiag, NULL, NULL, GrB_MINV_FP32, indiag, NULL));
        
        OK (GrB_Vector_new(&x, GrB_FP32, n));
        OK (GrB_assign(x, NULL, NULL, 1, GrB_ALL, n, NULL));
        OK (GrB_Vector_setElement_FP32(x, 0, 0));

        OK (LAGraph_mypcg2 (&steper, &kk, Y, u, alpha, indiag, x, .000001, 50,
            msg)) ;

        //--------------------------------------------------------------------------
        // try the LAGraph_Hdip_Fiedler algorithm
        //--------------------------------------------------------------------------

        //Set kmax = [20,50]
        OK (GrB_Vector_new(&kmax, GrB_FP32, 2));
        OK (GrB_Vector_setElement_FP32(kmax, 20, 0));
        OK (GrB_Vector_setElement_FP32(kmax, 50, 1));

        OK (LAGraph_Hdip_Fiedler (&iters, &lambda_result,
            &fiedler_vector, Y, infnorm, kmax, 0.000001, 0.000001, msg)) ;

        //--------------------------------------------------------------------------
        // check the results
        //--------------------------------------------------------------------------

        float err = difference (fiedler_vector, files [k].fiedler) ;
        TEST_CHECK (err < 1e-4) ;
        err = fabs(lambda_result - files [k].lambda) ;
        TEST_CHECK (err < 1e-4) ;

        //--------------------------------------------------------------------------
        // print the results
        //--------------------------------------------------------------------------

        printf("\n===============================The result vector x:\n");
        LAGraph_Vector_Print (fiedler_vector, 3, stdout, msg) ;
        printf("\n===============================The lambda: %f\n", lambda_result);
        printf("\n===============================The iters: \n");
        LAGraph_Vector_Print (iters, 3, stdout, msg) ;

        //--------------------------------------------------------------------------
        // free everyting and finish
        //--------------------------------------------------------------------------

        GrB_free (&A) ;
        GrB_free (&Y) ;
        GrB_free (&x) ;
        GrB_free (&fiedler_vector) ;
        GrB_free (&u) ;
        GrB_free (&steper) ;
        GrB_free (&indiag) ;
        GrB_free (&iters) ;
        GrB_free (&kmax) ;
    }

    OK (LAGraph_Finalize (msg)) ;
#endif
}

TEST_LIST = {
    {"Fieder", test_fiedler},
    {NULL, NULL}
} ;

