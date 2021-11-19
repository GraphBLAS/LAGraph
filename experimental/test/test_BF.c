//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_BF
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Contributed by Jinhao Chen and Tim Davis, Texas A&M

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>

//------------------------------------------------------------------------------
// globals
//------------------------------------------------------------------------------

#define LEN 512
char filename [LEN+1] ;
char msg [LAGRAPH_MSG_LEN] ;

//------------------------------------------------------------------------------
// test cases
//------------------------------------------------------------------------------

const char *files [ ] =
{
    "karate.mtx",
    "west0067.mtx",
    ""
} ;

//------------------------------------------------------------------------------
// setup: start a test
//------------------------------------------------------------------------------

void setup (void)
{
    OK (LAGraph_Init (msg)) ;
    OK (LAGraph_Random_Init (msg)) ;
}

//------------------------------------------------------------------------------
// teardown: finalize a test
//------------------------------------------------------------------------------

void teardown (void)
{
    OK (LAGraph_Random_Finalize (msg)) ;
    OK (LAGraph_Finalize (msg)) ;
}

#define LAGraph_FREE_ALL            \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&AT) ;                \
    free(I);                        \
    free(J);                        \
    free(W);                        \
    free(d);                        \
    free(pi);                       \
    GrB_free (&d1) ;                \
    GrB_free (&pi1) ;               \
    GrB_free (&h1) ;                \
    GrB_free (&d2) ;                \
    GrB_free (&pi2) ;               \
    GrB_free (&h2) ;                \
    GrB_free (&d3) ;                \
    GrB_free (&d4) ;                \
    GrB_free (&d5) ;                \
    GrB_free (&pi5) ;               \
    GrB_free (&h5) ;                \
    GrB_free (&d5a) ;               \
    GrB_free (&pi5a) ;              \
    GrB_free (&h5a) ;               \
    GrB_free (&d6) ;                \
    GrB_free (&d7) ;                \
    GrB_free (&pi6) ;               \
    GrB_free (&h6) ;                \
}

// int main (int argc, char **argv)

void test_BF (void)
{

GrB_Info info ;
setup ( ) ;

for (int k = 0 ; ; k++)
{

    double tic [2] ;
    GrB_Type A_type = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix AT = NULL ;
    GrB_Index *I = NULL, *J = NULL; // for col/row indices of entries from A
    double *W = NULL;
    double *d = NULL;
    int64_t *pi = NULL;
    GrB_Vector d1 = NULL ;
    GrB_Vector pi1 = NULL ;
    GrB_Vector h1 = NULL ;
    GrB_Vector d2 = NULL ;
    GrB_Vector pi2 = NULL ;
    GrB_Vector h2 = NULL ;
    GrB_Vector d3 = NULL ;
    GrB_Vector d4 = NULL ;
    GrB_Vector d5 = NULL ;
    GrB_Vector pi5 = NULL ;
    GrB_Vector h5 = NULL ;
    GrB_Vector d5a = NULL ;
    GrB_Vector pi5a = NULL ;
    GrB_Vector h5a = NULL ;
    GrB_Vector d6 = NULL ;
    GrB_Vector d7 = NULL ;
    GrB_Vector pi6 = NULL ;
    GrB_Vector h6 = NULL ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    // load the matrix as A
    const char *aname = files [k] ;
    if (strlen (aname) == 0) break;
    TEST_CASE (aname) ;
    snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &A_type, f, msg)) ;
    OK (fclose (f)) ;
    TEST_MSG ("Loading of valued matrix failed") ;
    printf ("\nMatrix: %s\n", aname) ;
    OK (LAGraph_Matrix_print (A, 2, stdout, NULL)) ;

    // TODO: all for some negative weights, and check for negative-weight cycles

    // ensure the matrix has all positive weights
    OK (GrB_apply (A, NULL, NULL, GrB_ABS_FP64, A, NULL)) ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nvals ;
    GrB_Matrix_nvals (&nvals, A) ;
    GrB_Index nrows, ncols ;
    OK (GrB_Matrix_nrows (&nrows, A)) ;
    OK (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;

    I = LAGraph_Malloc (nvals, sizeof(GrB_Index)) ;
    J = LAGraph_Malloc (nvals, sizeof(GrB_Index)) ;
    W = LAGraph_Malloc (nvals, sizeof(double)) ;

    OK (GrB_Matrix_extractTuples_FP64(I, J, W, &nvals, A));

    //--------------------------------------------------------------------------
    // get the source node
    //--------------------------------------------------------------------------

    GrB_Index s = 0 ;

    fprintf (stderr, "\n=========="
        "input graph: nodes: %lu edges: %lu source node: %lu\n", n, nvals, s) ;

    //--------------------------------------------------------------------------
    // run LAGraph_BF_full1 before setting the diagonal to 0
    //--------------------------------------------------------------------------

    int ntrials = 1 ;       // increase this to 10, 100, whatever, for more
                            // accurate timing
    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d5) ;
        GrB_free (&pi5) ;
        GrB_free (&h5) ;
        int result = (LAGraph_BF_full1 (&d5, &pi5, &h5, A, s)) ;
        printf ("result: %d\n", result) ;
        OK (result) ;
    }

    // stop the timer
    double t5;
    OK (LAGraph_Toc (&t5, tic, NULL)) ;
    t5 = t5 / ntrials;
    fprintf (stderr, "BF_full1      time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t5, 1e-6*((double) nvals) / t5) ;
    //--------------------------------------------------------------------------
    // run LAGraph_BF_full1a before setting the diagonal to 0
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d5a) ;
        GrB_free (&pi5a) ;
        GrB_free (&h5a) ;
        OK (LAGraph_BF_full1a (&d5a, &pi5a, &h5a, A, s)) ;
    }

    // stop the timer
    double t5a;
    OK (LAGraph_Toc (&t5a, tic, NULL)) ;
    t5a = t5a / ntrials;
    fprintf (stderr, "BF_full1a     time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t5a, 1e-6*((double) nvals) / t5a) ;

    //--------------------------------------------------------------------------
    // run LAGraph_BF_full2 before setting the diagonal to 0
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d6) ;
        GrB_free (&pi6) ;
        GrB_free (&h6) ;
        OK (LAGraph_BF_full2 (&d6, &pi6, &h6, A, s)) ;
    }

    // stop the timer
    double t6;
    OK (LAGraph_Toc (&t6, tic, NULL)) ;
    t6 = t6 / ntrials;
    fprintf (stderr, "BF_full2      time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t6, 1e-6*((double) nvals) / t6) ;


    //--------------------------------------------------------------------------
    // set the diagonal to 0
    //--------------------------------------------------------------------------
    for (GrB_Index i = 0; i < n; i++)
    {
        OK (GrB_Matrix_setElement_FP64 (A, 0, i, i));
    }

    //--------------------------------------------------------------------------
    // AT = A'
    //--------------------------------------------------------------------------

    OK (LAGraph_Tic (tic, NULL)) ;
    OK (GrB_Matrix_new (&AT, GrB_FP64, ncols, nrows)) ;
    OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
    double transpose_time;
    OK (LAGraph_Toc (&transpose_time, tic, NULL)) ;
    fprintf (stderr, "transpose     time: %g\n", transpose_time) ;

    //--------------------------------------------------------------------------
    // run the LAGraph_BF_full on node s
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d1) ;
        GrB_free (&pi1) ;
        GrB_free (&h1) ;
        OK (LAGraph_BF_full (&d1, &pi1, &h1, A, s)) ;
    }

    // stop the timer
    double t1;
    OK (LAGraph_Toc (&t1, tic, NULL)) ;
    t1 = t1 / ntrials;
    fprintf (stderr, "BF_full       time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t1, 1e-6*((double) nvals) / t1) ;
    fprintf (stderr, "t(BF_full1) / t(BF_full):      %g\n", t5/t1) ;

    //--------------------------------------------------------------------------
    // run the BF on node s with LAGraph_BF_basic
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d3) ;
        OK (LAGraph_BF_basic (&d3, A, s)) ;
    }

    // stop the timer
    double t2;
    OK (LAGraph_Toc (&t2, tic, NULL)) ;
    t2 = t2 / ntrials;
    fprintf (stderr, "BF_basic      time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t2, 1e-6*((double) nvals) / t2) ;
    fprintf (stderr, "speedup of BF_basic:       %g\n", t1/t2) ;

    //--------------------------------------------------------------------------
    // run the BF on node s with LAGraph_pure_c
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        free (d) ;
        free (pi) ;
        OK (LAGraph_BF_pure_c_double (&d, &pi, s, n, nvals, I, J, W)) ;
    }

    // stop the timer
    double t3;
    OK (LAGraph_Toc (&t3, tic, NULL)) ;
    t3 = t3 / ntrials;
    fprintf (stderr, "BF_pure_c     time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t3, 1e-6*((double) nvals) / t3) ;
    fprintf (stderr, "speedup of BF_pure_c:      %g\n", t1/t3) ;

    //--------------------------------------------------------------------------
    // run the LAGraph_BF_full_mxv on node s
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d2) ;
        GrB_free (&pi2) ;
        GrB_free (&h2) ;
        OK (LAGraph_BF_full_mxv (&d2, &pi2, &h2, AT, s)) ;
    }

    // stop the timer
    double t4;
    OK (LAGraph_Toc (&t4, tic, NULL)) ;
    t4 = t4 / ntrials;
    fprintf (stderr, "BF_full_mxv   time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t4, 1e-6*((double) nvals) / t4) ;
    fprintf (stderr, "speedup of BF_full_mxv:    %g\n", t1/t4) ;

    //--------------------------------------------------------------------------
    // run the BF on node s with LAGraph_BF_basic_mxv
    //--------------------------------------------------------------------------

    // start the timer
    OK (LAGraph_Tic (tic, NULL)) ;

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_free (&d4) ;
        OK (LAGraph_BF_basic_mxv (&d4, AT, s)) ;
    }

    // stop the timer
    double t;
    OK (LAGraph_Toc (&t5, tic, NULL)) ;
    t5 = t5 / ntrials;
    fprintf (stderr, "BF_basic_mxv  time: %12.6e (sec), rate:"
        " %g (1e6 edges/sec)\n", t5, 1e-6*((double) nvals) / t5) ;
    fprintf (stderr, "speedup of BF_basic_mxv:   %g\n", t1/t5) ;

    //--------------------------------------------------------------------------
    // run the BF on node s with LAGraph_BF_basic_pushpull
    //--------------------------------------------------------------------------

    GrB_free (&d7) ;
    OK (LAGraph_BF_basic_pushpull (&d7, A, AT, s)) ;

    //--------------------------------------------------------------------------
    // check results
    //--------------------------------------------------------------------------

    bool isequal = false ;

    TEST_CHECK (d != NULL && d1 != NULL) ;
    {
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = INFINITY ;
            int64_t pii = 0;
            OK (GrB_Vector_extractElement (&di, d1, i)) ;
            TEST_CHECK (di == d[i]) ;

            // since d5 is a dense vector filled with infinity, we have to
            // compare it against d seperaterly
            OK (GrB_Vector_extractElement (&di, d5, i)) ;
            TEST_CHECK (di == d[i]) ;

            // since d5a is a dense vector filled with infinity, we have to
            // compare it against d seperaterly
            OK (GrB_Vector_extractElement (&di, d5a, i)) ;
            TEST_CHECK (di == d[i]) ;
            /*
            OK (GrB_Vector_extractElement (&pii, pi1, i)) ;
            if (pii != pi[i]+1)
            {
                fprintf (stderr, "ERROR! BF_full and BF_pure_c pi  differ\n") ;
                TEST_CHECK (false) ;
                break;
            }
            */
        }
    }
//  else
//  {
//      fprintf (stderr, "ERROR! BF_full and BF_pure_c d  differ\n") ;
//      fprintf (stderr, "BF_full %s negative-weight cycle, while BF_pure_c %s "
//          "negative-weight cycle\n", (d1 == NULL) ? "found":"didn't find",
//          (d == NULL) ? "found":"didn't find");
//      TEST_CHECK (false) ;
//  }

    OK (LAGraph_Vector_IsEqual (&isequal, d1, d3, NULL)) ;
    TEST_CHECK (isequal) ;
//  if (!isequal)
//  {
//      fprintf (stderr, "ERROR! BF_full and BF_basic   differ\n") ;
//      TEST_CHECK (false) ;
//  }
    OK (LAGraph_Vector_IsEqual (&isequal, d1, d4, NULL)) ;
    TEST_CHECK (isequal) ;
//  if (!isequal)
//  {
//      fprintf (stderr, "ERROR! BF_full and BF_basic_mxv   differ\n") ;
//      TEST_CHECK (false) ;
//  }
    OK (LAGraph_Vector_IsEqual (&isequal, d1, d2, NULL)) ;
    TEST_CHECK (isequal) ;
//  if (!isequal)
//  {
//      fprintf (stderr, "ERROR! BF_full and BF_full_mxv    differ\n") ;
//      TEST_CHECK (false) ;
//  }

    OK (LAGraph_Vector_IsEqual (&isequal, d1, d6, NULL)) ;
    TEST_CHECK (isequal) ;
//  if (!isequal)
//  {
//      fprintf (stderr, "ERROR! BF_full and BF_full2   differ\n") ;
//      TEST_CHECK (false) ;
//  }
/*
    OK (LAGraph_Vector_IsEqual (&isequal, pi1, pi2, NULL)) ;
    if (!isequal)
    {
        fprintf (stderr, "ERROR! BF_full and BF_full_mxv pi   differ\n") ;
        TEST_CHECK (false) ;
    }
*/

    OK (LAGraph_Vector_IsEqual (&isequal, d1, d7, NULL)) ;
    TEST_CHECK (isequal) ;

    //--------------------------------------------------------------------------
    // write the result to stdout (check them outside of this main program)
    //--------------------------------------------------------------------------

#if 0
    fprintf (stderr,
    "\n------------------------------------------------------------\n\n") ;
    // pure c
    if (d != NULL)
    {
        fprintf (stderr, "BF pure C result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            fprintf (stderr, "%4.2f  ", d[i]) ;
        }
        fprintf (stderr, "\n") ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            fprintf (stderr, "%4" PRIu64 "  ", pi[i]+1) ;
        }
        fprintf (stderr, "\n") ;
        fprintf (stderr, "\n") ;
    }
    // BF full1
    if (d5 != NULL)
    {
        fprintf (stderr, "BF full1 result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = 0 ;
            OK (GrB_Vector_extractElement (&di, d5, i)) ;
            fprintf (stderr, "%4.2f  ", di) ;
        }
        fprintf (stderr, "\n") ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            int64_t x = 0 ;
            OK (GrB_Vector_extractElement (&x, pi5, i)) ;
            fprintf (stderr, "%4" PRIu64 "  ", x) ;
        }
        fprintf (stderr, "\n") ;
        fprintf (stderr, "\n") ;
    }
    // BF full1a
    if (d5a != NULL)
    {
        fprintf (stderr, "BF full1a result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = 0 ;
            OK (GrB_Vector_extractElement (&di, d5a, i)) ;
            fprintf (stderr, "%4.2f  ", di) ;
        }
        fprintf (stderr, "\n") ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            int64_t x = 0 ;
            OK (GrB_Vector_extractElement (&x, pi5a, i)) ;
            fprintf (stderr, "%4" PRIu64 "  ", x) ;
        }
        fprintf (stderr, "\n") ;
        fprintf (stderr, "\n") ;
    }
    // BF full
    if (d1 != NULL)
    {
        fprintf (stderr, "BF full result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = 0 ;
            OK (GrB_Vector_extractElement (&di, d1, i)) ;
            fprintf (stderr, "%4.2f  ", di) ;
        }
        fprintf (stderr, "\n") ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            int64_t x = 0 ;
            OK (GrB_Vector_extractElement (&x, pi1, i)) ;
            fprintf (stderr, "%4" PRIu64 "  ", x) ;
        }
        fprintf (stderr, "\n") ;
        fprintf (stderr, "\n") ;
    }
    // BF full mxv
    if (d2 != NULL)
    {
        fprintf (stderr, "BF full mxv result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = 0 ;
            OK (GrB_Vector_extractElement (&di, d2, i)) ;
            fprintf (stderr, "%4.2f  ", di) ;
        }
        printf ("\n") ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            int64_t x = 0 ;
            OK (GrB_Vector_extractElement (&x, pi2, i)) ;
            fprintf (stderr, "%4" PRIu64 "  ", x) ;
        }
        printf ("\n") ;
        printf ("\n") ;
    }
    // BF basic
    if (d3 != NULL)
    {
        fprintf (stderr, "BF basic result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = 0 ;
            OK (GrB_Vector_extractElement (&di, d3, i)) ;
            fprintf (stderr, "%4.2f  ", di) ;
        }
        fprintf (stderr, "\n") ;
        fprintf (stderr, "\n") ;
    }
    if (d4 != NULL)
    {
        fprintf (stderr, "BF basic mxv result\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            double di = 0 ;
            OK (GrB_Vector_extractElement (&di, d4, i)) ;
            fprintf (stderr, "%4.2f  ", di) ;
        }
        fprintf (stderr, "\n") ;
    }

    if (d != NULL)
    {
        fprintf (stderr, "shortest path for each node\n");
        for (int64_t i = 0 ; i < n ; i++)
        {
            int64_t node = i ;
            fprintf (stderr, "%4" PRIu64 "", node+1) ;
            while (node != s && node != -1 )
            {
                //OK (GrB_Vector_extractElement (&node, pi, node-1)) ;
                //printf (" <- %" PRIu64 "", node) ;
                node = pi[node] ;
                fprintf (stderr, " <- %4" PRIu64 "", node+1) ;
            }
            fprintf (stderr, "\n") ;
        }
    }
#endif

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGraph_FREE_ALL ;
}

teardown ( ) ;
}

//------------------------------------------------------------------------------
// TEST_LIST: list of tasks for this entire test
//------------------------------------------------------------------------------

TEST_LIST =
{
    { "test_BF", test_BF },
    { NULL, NULL }
} ;

