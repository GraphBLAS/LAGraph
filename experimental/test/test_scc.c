//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_scc.c: tests for Strongly Connected Components
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// todo: write a simple scc method, as LG_check_scc, and compare its results
// with LAGraph_scc

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    const char *name ;
    int cc_count;
    uint64_t hash;
}
matrix_info ;

int scc_cover [7] = { 0, 0, 2, 0, 4, 2, 0 } ;

const matrix_info files [ ] =
{
    { "A2.mtx", 1, 0x2de8d717be626313},
    { "A.mtx", 1, 0x2de8d717be626313},
    { "bcsstk13.mtx", 1, 0x41d903f08b46b543},
    { "cover.mtx", 3, 0x30ae8cb78a807691},
    { "cover_structure.mtx", 3, 0x30ae8cb78a807691},
    { "cryg2500.mtx", 1, 0xd1cb8e3cc6be967},
    { "full.mtx", 1, 0x99971e4f016b4644},
    { "full_noheader.mtx", 1, 0x99971e4f016b4644},
    { "full_symmetric.mtx", 1, 0x278859fec1de1f7f},
    { "jagmesh7.mtx", 1, 0x66b315eea17941c8},
    { "karate.mtx", 1, 0x8bad7c50644c4aa9},
    { "ldbc-cdlp-directed-example.mtx", 2, 0x3a61ac294b7bb114},
    { "ldbc-cdlp-undirected-example.mtx", 1, 0x4072e255fd8e310a},
    { "ldbc-directed-example-bool.mtx", 7, 0xc66f5ecf1b7f6876},
    { "ldbc-directed-example.mtx", 7, 0xc66f5ecf1b7f6876},
    { "ldbc-directed-example-unweighted.mtx", 7, 0xc66f5ecf1b7f6876},
    { "ldbc-undirected-example-bool.mtx", 1, 0xf53db7dbbeff3283},
    { "ldbc-undirected-example.mtx", 1, 0xf53db7dbbeff3283},
    { "ldbc-undirected-example-unweighted.mtx", 1, 0xf53db7dbbeff3283},
    { "ldbc-wcc-example.mtx", 1, 0x36a78022528a2101},
    { "LFAT5.mtx", 3, 0x79d4d8de0a22a863},
    { "LFAT5_two.mtx", 6, 0xac369d0362d73d6},
    { "matrix_bool.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_fp32.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_fp32_structure.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_fp64.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_int16.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_int32.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_int64.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_int8.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_uint16.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_uint32.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_uint64.mtx", 3, 0x30ae8cb78a807691},
    { "matrix_uint8.mtx", 3, 0x30ae8cb78a807691},
    { "msf1.mtx", 4, 0x6445d984131a9555},
    { "msf2.mtx", 8, 0x72d532720c54b673},
    { "msf3.mtx", 5, 0xf57eb057beb5a5c7},
    { "olm1000.mtx", 1, 0x15cf9ea2db88ab18},
    { "pushpull.mtx", 1, 0x1816384cd04f7e01},
    { "sample2.mtx", 1, 0x4072e255fd8e310a},
    { "sample.mtx", 8, 0x72d532720c54b673},
    { "structure.mtx", 3, 0x30ae8cb78a807691},
    { "test_BF.mtx", 3, 0x30ae8cb78a807691},
    { "test_FW_1000.mtx", 1, 0x15cf9ea2db88ab18},
    { "test_FW_2003.mtx", 485, 0xf79ad45d3a704eec},
    { "test_FW_2500.mtx", 646, 0x4fa83d60352e7e19},
    { "tree-example.mtx", 1, 0x8857b82baeba129},
    { "west0067_jumbled.mtx", 1, 0xa861dc7526128ac7},
    { "west0067.mtx", 1, 0xa861dc7526128ac7},
    { "west0067_noheader.mtx", 1, 0xa861dc7526128ac7},
    { "zenios.mtx", 1391, 0x15b2b99a80c3480e},
    { "", 0, 0},
} ;

//------------------------------------------------------------------------------
// count_connected_components: count the # of components in a component vector
//------------------------------------------------------------------------------

int count_connected_components (GrB_Vector C)
{
    GrB_Index n = 0 ;
    OK (GrB_Vector_size (&n, C)) ;
    int ncomponents = 0 ;
    for (int i = 0 ; i < n ; i++)
    {
        int64_t comp = -1 ;
        int result = GrB_Vector_extractElement (&comp, C, i) ;
        if (result == GrB_SUCCESS && comp == i) ncomponents++ ;
    }
    return (ncomponents) ;
}

//****************************************************************************
void test_scc (void)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LAGraph_Init (msg) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;

        GrB_Vector c = NULL ;

        for (int jit = 0 ; jit <= 1 ; jit++)
        {
            OK (GxB_Global_Option_set (GxB_JIT_C_CONTROL,
                jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;
            // find the strongly connected components with LAGraph_scc
            // GrB_set (GrB_GLOBAL, (int32_t) (true), GxB_BURBLE) ;
            OK (LAGraph_scc (&c, A, msg)) ;
            // GrB_set (GrB_GLOBAL, (int32_t) (true), GxB_BURBLE) ;

            GrB_Index n ;
            OK (GrB_Vector_size (&n, c)) ;
            LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT ;

            // check result c for cover
            if (strcmp (aname, "cover.mtx") == 0)
            {
                GrB_Vector cgood = NULL ;
                OK (GrB_Vector_new (&cgood, GrB_UINT64, n)) ;
                for (int k = 0 ; k < n ; k++)
                {
                    OK (GrB_Vector_setElement (cgood, scc_cover [k], k)) ;
                }
                OK (GrB_wait (cgood, GrB_MATERIALIZE)) ;
                printf ("\nscc (known result):\n") ;
                OK (LAGraph_Vector_Print (cgood, pr, stdout, msg)) ;
                bool ok = false ;
                OK (LAGraph_Vector_IsEqual (&ok, c, cgood, msg)) ;
                TEST_CHECK (ok) ;
                OK (GrB_free (&cgood)) ;
            }
            int result_cc_count = count_connected_components(c);
            TEST_CHECK(result_cc_count == files[k].cc_count);
            OK (LAGraph_Vector_Print (c, pr, stdout, msg)) ;
            uint64_t hash = 0;
            int hash_info = LAGraph_Hash_Vector(&hash, c, msg); 
            OK (hash_info);
            TEST_CHECK(hash == files[k].hash);
            OK (GrB_free (&c)) ;
        }
        OK (GrB_free (&A)) ;
    }

    LAGraph_Finalize (msg) ;
    #endif
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors (void)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    LAGraph_Init (msg) ;

    GrB_Vector c = NULL ;
    GrB_Matrix A = NULL ;

    // c and A are NULL
    int result = LAGraph_scc (NULL, A, msg) ;
    printf ("\nresult: %d\n", result) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // A is rectangular
    OK (GrB_Matrix_new (&A, GrB_BOOL, 3, 4)) ;
    result = LAGraph_scc (&c, A, msg) ;
    TEST_CHECK (result == GrB_DIMENSION_MISMATCH) ;

    OK (GrB_free (&c)) ;
    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
    #endif
}

//****************************************************************************

TEST_LIST = {
    {"scc", test_scc},
    {"scc_errors", test_errors},
    {NULL, NULL}
};
