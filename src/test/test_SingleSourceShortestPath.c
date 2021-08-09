//----------------------------------------------------------------------------
// LAGraph/src/test/test_SingleSourceShortestPath.cpp: test cases for triangle
// counting algorithms
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
        // GxB_print (A, 3) ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        OK (LAGraph_CheckGraph (G, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // delta values to try
        int32_t Deltas [ ] = { 30, 100, 50000 } ;

        // run the SSSP
        int64_t step = (n > 100) ? (3*n/4) : ((n/4) + 1) ;
        for (int64_t src = 0 ; src < n ; src += step)
        {
            GrB_Vector path_length = NULL ;
            for (int32_t kk = 0 ; kk < 3 ; kk++)
            {
                int32_t delta = Deltas [kk] ;
                printf ("\n========================= src %ld delta %d n %ld\n",
                    src, delta, n) ;
                OK (LAGraph_SingleSourceShortestPath (&path_length,
                    G, src, delta, true, msg)) ;
                int result = LG_check_sssp (path_length, G, src, msg) ;
                if (result != 0)
                {
                    printf ("failure %d %s\n", result, msg) ;
                }
                OK (result) ;
                OK (GrB_free(&path_length)) ;
            }
        }
        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize(msg);
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"SSSP", test_SingleSourceShortestPath},
    {NULL, NULL}
};

