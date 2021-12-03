//----------------------------------------------------------------------------
// LAGraph/src/test/test_SingleSourceShortestPath.cpp: test cases for SSSP
// ----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;

#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] = 
{
    { "A.mtx" }, 
    { "cover.mtx" }, 
    { "jagmesh7.mtx" }, 
    { "ldbc-cdlp-directed-example.mtx" }, 
    { "ldbc-cdlp-undirected-example.mtx" }, 
    { "ldbc-directed-example.mtx" }, 
    { "ldbc-undirected-example.mtx" }, 
    { "ldbc-wcc-example.mtx" }, 
    { "LFAT5.mtx" }, 
    { "msf1.mtx" }, 
    { "msf2.mtx" }, 
    { "msf3.mtx" }, 
    { "sample2.mtx" }, 
    { "sample.mtx" }, 
    { "olm1000.mtx" }, 
    { "bcsstk13.mtx" }, 
    { "cryg2500.mtx" }, 
    { "tree-example.mtx" }, 
    { "west0067.mtx" }, 
    { "karate.mtx" }, 
    { "matrix_bool.mtx" }, 
    { "test_BF.mtx" },
    { "test_FW_1000.mtx" },
    { "test_FW_2003.mtx" },
    { "test_FW_2500.mtx" },
    { "skew_fp32.mtx" }, 
    { "" }, 
} ;

//****************************************************************************
void test_SingleSourceShortestPath(void)
{
    LAGraph_Init(msg);
    GrB_Matrix A = NULL, T = NULL ;
    GrB_Type atype = NULL ;

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        LAGraph_Kind kind = LAGRAPH_ADJACENCY_DIRECTED ;
        // LAGraph_Kind kind = files [k].kind ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\nMatrix: %s\n", aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        GrB_Index n = 0, ncols = 1 ;
        OK (GrB_Matrix_nrows (&n, A)) ;
        OK (GrB_Matrix_ncols (&ncols, A)) ;
        TEST_CHECK (n == ncols) ;

        // convert A to int32
        if (atype != GrB_INT32)
        {
            OK (GrB_Matrix_new (&T, GrB_INT32, n, n)) ;
            OK (GrB_assign (T, NULL, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL)) ;
            atype = GrB_INT32 ;
            OK (GrB_free (&A)) ;
            A = T ;
        }

        // ensure all entries are positive, and in the range 1 to 255
        OK (GrB_Matrix_apply_BinaryOp2nd_INT32 (A, NULL, NULL,
            GrB_BAND_INT32, A, 255, NULL)) ;
        int32_t x ; 
        OK (GrB_reduce (&x, NULL, GrB_MIN_MONOID_INT32, A, NULL)) ;
        if (x < 1)
        {
            OK (GrB_Matrix_apply_BinaryOp2nd_INT32 (A, NULL, NULL,
                GrB_MAX_INT32, A, 1, NULL)) ;
        }

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // delta values to try
        int32_t Deltas [ ] = { 30, 100, 50000 } ;

        // run the SSSP
        GrB_Vector path_length = NULL ;
        int64_t step = (n > 100) ? (3*n/4) : ((n/4) + 1) ;
        for (int64_t src = 0 ; src < n ; src += step)
        {
            for (int32_t kk = 0 ; kk < ((n > 100) ? 1 : 3) ; kk++)
            {
                int32_t delta = Deltas [kk] ;
                printf ("src %d delta %d n %d\n", (int) src, delta, (int) n) ;
                OK (LAGraph_SingleSourceShortestPath (&path_length,
                    G, src, delta, true, msg)) ;
                OK (LG_check_sssp (path_length, G, src, msg)) ;
                OK (GrB_free(&path_length)) ;
            }
        }

        // add a single negative edge and try again
        OK (GrB_Matrix_setElement_INT32 (G->A, -1, 0, 1)) ;
        OK (LAGraph_SingleSourceShortestPath (&path_length,
            G, 0, 30, false, msg)) ;
        OK (LAGraph_Vector_print (path_length, 2, stdout, msg)) ;
        int32_t len = 0 ;
        OK (GrB_Vector_extractElement (&len, path_length, 1)) ;
        TEST_CHECK (len == -1) ;
        OK (GrB_free(&path_length)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize(msg);
}

//------------------------------------------------------------------------------
// test_SingleSourceShortestPath_brutal
//------------------------------------------------------------------------------

#if LG_SUITESPARSE
void test_SingleSourceShortestPath_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;

    GrB_Matrix A = NULL, T = NULL ;
    GrB_Type atype = NULL ;

    // just test with the first 8 matrices
    for (int k = 0 ; k < 8 ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        LAGraph_Kind kind = LAGRAPH_ADJACENCY_DIRECTED ;
        // LAGraph_Kind kind = files [k].kind ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\nMatrix: %s\n", aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        GrB_Index n = 0 ;
        OK (GrB_Matrix_nrows (&n, A)) ;
        if (n > 30)
        {
            printf ("skipped -- only using small matrices for brutal test\n") ;
            OK (GrB_free (&A)) ;
            continue ;
        }

        // convert A to int32
        if (atype != GrB_INT32)
        {
            OK (GrB_Matrix_new (&T, GrB_INT32, n, n)) ;
            OK (GrB_assign (T, NULL, NULL, A, GrB_ALL, n, GrB_ALL, n, NULL)) ;
            atype = GrB_INT32 ;
            OK (GrB_free (&A)) ;
            A = T ;
        }

        // ensure all entries are positive, and in the range 1 to 255
        OK (GrB_Matrix_apply_BinaryOp2nd_INT32 (A, NULL, NULL,
            GrB_BAND_INT32, A, 255, NULL)) ;
        int32_t x ; 
        OK (GrB_reduce (&x, NULL, GrB_MIN_MONOID_INT32, A, NULL)) ;
        if (x < 1)
        {
            OK (GrB_Matrix_apply_BinaryOp2nd_INT32 (A, NULL, NULL,
                GrB_MAX_INT32, A, 1, NULL)) ;
        }

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;

        // run the SSSP on a single source node with one delta
        GrB_Vector path_length = NULL ;
        int64_t src = 0 ;
        int32_t delta = 30 ;
        printf ("src %d delta %d n %d\n", (int) src, delta, (int) n) ;
        LG_BRUTAL (LAGraph_SingleSourceShortestPath (&path_length,
            G, src, delta, true, msg)) ;
        OK (LG_check_sssp (path_length, G, src, msg)) ;
        OK (GrB_free(&path_length)) ;

        // add a single negative edge and try again
        OK (GrB_Matrix_setElement_INT32 (G->A, -1, 0, 1)) ;
        OK (GrB_wait (G->A, GrB_MATERIALIZE)) ;
        LG_BRUTAL (LAGraph_SingleSourceShortestPath (&path_length,
            G, 0, 30, false, msg)) ;
        int32_t len = 0 ;
        OK (GrB_Vector_extractElement (&len, path_length, 1)) ;
        TEST_CHECK (len == -1) ;
        OK (GrB_free(&path_length)) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//****************************************************************************
//****************************************************************************

TEST_LIST = {
    {"SSSP", test_SingleSourceShortestPath},
    #if LG_SUITESPARSE
    {"SSSP_brutal", test_SingleSourceShortestPath_brutal },
    #endif
    {NULL, NULL}
};

