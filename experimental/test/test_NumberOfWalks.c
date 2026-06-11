//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_NumberOfWalks.c: test for NumberOfWalks
//------------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>
#include "LG_Xtest.h"

#define LEN 512
char msg [LAGRAPH_MSG_LEN] ;
char filename [LEN+1] ;

//------------------------------------------------------------------------------
// Ground truth results (pre-computed with NetworkX)
//------------------------------------------------------------------------------

// A.mtx (7x7), A^2 walks of length 2
int64_t A_mtx_walks_k2[49] = {
    3, 1, 2, 3, 2, 3, 1,
    1, 5, 3, 3, 1, 2, 5,
    2, 3, 5, 3, 3, 3, 3,
    3, 3, 3, 5, 2, 3, 3,
    2, 1, 3, 2, 3, 3, 1,
    3, 2, 3, 3, 3, 4, 2,
    1, 5, 3, 3, 1, 2, 5
} ;

// 4-node cycle, A^2 walks of length 2
int64_t cycle_4node_walks_k2[16] = {
    2, 0, 2, 0,
    0, 2, 0, 2,
    2, 0, 2, 0,
    0, 2, 0, 2
} ;

// 3-node path (0-1-2), A^2 walks of length 2
int64_t path_3node_walks_k2[9] = {
    1, 0, 1,
    0, 2, 0,
    1, 0, 1
} ;

// check matrix against ground truth
int64_t check_walks (GrB_Matrix result, int64_t *truth, GrB_Index n)
{
    int64_t max_diff = 0 ;
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        for (GrB_Index j = 0 ; j < n ; j++)
        {
            int64_t computed = 0 ;
            int64_t expected = truth [i * n + j] ;
            
            GrB_Info info = GrB_Matrix_extractElement_INT64 (&computed, result, i, j) ;
            if (info == GrB_NO_VALUE) computed = 0 ;
            
            int64_t diff = computed > expected ? computed - expected : expected - computed ;
            if (diff > max_diff) max_diff = diff ;
        }
    }
    return max_diff ;
}

//------------------------------------------------------------------------------
// test_NumberOfWalks_path_3node: simple path graph
//------------------------------------------------------------------------------

void test_NumberOfWalks_path_3node (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, result = NULL ;
    GrB_Index n = 3 ;
    
    OK (GrB_Matrix_new (&A, GrB_INT64, n, n)) ;
    
    GrB_Index rows[] = {0, 1, 1, 2} ;
    GrB_Index cols[] = {1, 0, 2, 1} ;
    int64_t vals[] = {1, 1, 1, 1} ;
    
    OK (GrB_Matrix_build_INT64 (A, rows, cols, vals, 4, GrB_PLUS_INT64)) ;
    
    // Compute A^2
    OK (LAGraph_NumberOfWalks (&result, A, NULL, 2)) ;

    // Compare with ground truth
    int64_t err = check_walks (result, path_3node_walks_k2, n) ;
    TEST_CHECK (err == 0) ;
    
    OK (GrB_free (&A)) ;
    OK (GrB_free (&result)) ;
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_NumberOfWalks_A_mtx: A.mtx graph with k=2
//------------------------------------------------------------------------------

void test_NumberOfWalks_A_mtx (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, result = NULL ;
    GrB_Index n = 7 ;
    
    snprintf (filename, LEN, LG_DATA_DIR "%s", "A.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    
    OK (LAGraph_NumberOfWalks (&result, A, NULL, 2)) ;

    int64_t err = check_walks (result, A_mtx_walks_k2, n) ;
    TEST_CHECK (err == 0) ;
    
    OK (GrB_free (&A)) ;
    OK (GrB_free (&result)) ;
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_NumberOfWalks_cycle_4: simple 4-cycle
//------------------------------------------------------------------------------

void test_NumberOfWalks_cycle_4 (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, result = NULL ;
    GrB_Index n = 4 ;
    
    // Build 4-cycle: 0-1-2-3-0
    OK (GrB_Matrix_new (&A, GrB_INT64, n, n)) ;
    
    GrB_Index rows[] = {0, 1, 1, 2, 2, 3, 3, 0} ;
    GrB_Index cols[] = {1, 0, 2, 1, 3, 2, 0, 3} ;
    int64_t vals[] = {1, 1, 1, 1, 1, 1, 1, 1} ;
    
    OK (GrB_Matrix_build_INT64 (A, rows, cols, vals, 8, GrB_PLUS_INT64)) ;
    OK (LAGraph_NumberOfWalks (&result, A, NULL, 2)) ;
    
    int64_t err = check_walks (result, cycle_4node_walks_k2, n) ;
    TEST_CHECK (err == 0) ;
    
    OK (GrB_free (&A)) ;
    OK (GrB_free (&result)) ;
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_NumberOfWalks_varying_k: different walk lengths
//------------------------------------------------------------------------------

void test_NumberOfWalks_varying_k (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL, result = NULL ;
    GrB_Index n = 3 ;

    OK (GrB_Matrix_new (&A, GrB_INT64, n, n)) ;
    
    GrB_Index rows[] = {0, 1, 1, 2} ;
    GrB_Index cols[] = {1, 0, 2, 1} ;
    int64_t vals[] = {1, 1, 1, 1} ;
    
    OK (GrB_Matrix_build_INT64 (A, rows, cols, vals, 4, GrB_PLUS_INT64)) ;
    
    // k=1
    OK (LAGraph_NumberOfWalks (&result, A, NULL, 1)) ;
    int64_t val = 0 ;
    GrB_Matrix_extractElement_INT64 (&val, result, 0, 1) ;
    TEST_CHECK (val == 1) ;
    OK (GrB_free (&result)) ;
    
    // k=3
    OK (LAGraph_NumberOfWalks (&result, A, NULL, 3)) ;
    val = 0 ;
    GrB_Matrix_extractElement_INT64 (&val, result, 0, 2) ;
    TEST_CHECK (val == 0) ;
    OK (GrB_free (&result)) ;
    
    // k=4
    OK (LAGraph_NumberOfWalks (&result, A, NULL, 4)) ;
    val = 0 ;
    GrB_Matrix_extractElement_INT64 (&val, result, 0, 0) ;
    TEST_CHECK (val == 2) ;
    OK (GrB_free (&result)) ;
    
    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}


TEST_LIST = {
    {"test_NumberOfWalks_path_3node",  test_NumberOfWalks_path_3node},
    {"test_NumberOfWalks_A_mtx",       test_NumberOfWalks_A_mtx},
    {"test_NumberOfWalks_cycle_4",     test_NumberOfWalks_cycle_4},
    {"test_NumberOfWalks_varying_k",   test_NumberOfWalks_varying_k},
    {NULL, NULL}
} ;
