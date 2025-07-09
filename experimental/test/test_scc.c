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
    { "A2.mtx", 1, 6493938657738929428ull},
    { "A.mtx", 1, 6493938657738929428ull},
    { "bcsstk13.mtx", 1, 4873650117803742346ull},
    { "cover.mtx", 3, 848279640410529436ull},
    { "cover_structure.mtx", 3, 848279640410529436ull},
    { "cryg2500.mtx", 1, 8070599988610413093ull},
    { "full.mtx", 1, 15769435293242772098ull},
    { "full_noheader.mtx", 1, 15769435293242772098ull},
    { "full_symmetric.mtx", 1, 7920595475144714245ull},
    { "jagmesh7.mtx", 1, 5114200449021899176ull},
    { "karate.mtx", 1, 4176608668907330736ull},
    { "ldbc-cdlp-directed-example.mtx", 2, 11183292771650049706ull},
    { "ldbc-cdlp-undirected-example.mtx", 1, 4918287057298807835ull},
    { "ldbc-directed-example-bool.mtx", 7, 11304580677056001228ull},
    { "ldbc-directed-example.mtx", 7, 11304580677056001228ull},
    { "ldbc-directed-example-unweighted.mtx", 7, 11304580677056001228ull},
    { "ldbc-undirected-example-bool.mtx", 1, 9158223257798130275ull},
    { "ldbc-undirected-example.mtx", 1, 9158223257798130275ull},
    { "ldbc-undirected-example-unweighted.mtx", 1, 9158223257798130275ull},
    { "ldbc-wcc-example.mtx", 1, 4317729120311459500ull},
    { "LFAT5.mtx", 3, 17553140753101484131ull},
    { "LFAT5_two.mtx", 6, 7979561620824911ull},
    { "matrix_bool.mtx", 3, 848279640410529436ull},
    { "matrix_fp32.mtx", 3, 848279640410529436ull},
    { "matrix_fp32_structure.mtx", 3, 848279640410529436ull},
    { "matrix_fp64.mtx", 3, 848279640410529436ull},
    { "matrix_int16.mtx", 3, 848279640410529436ull},
    { "matrix_int32.mtx", 3, 848279640410529436ull},
    { "matrix_int64.mtx", 3, 848279640410529436ull},
    { "matrix_int8.mtx", 3, 848279640410529436ull},
    { "matrix_uint16.mtx", 3, 848279640410529436ull},
    { "matrix_uint32.mtx", 3, 848279640410529436ull},
    { "matrix_uint64.mtx", 3, 848279640410529436ull},
    { "matrix_uint8.mtx", 3, 848279640410529436ull},
    { "msf1.mtx", 4, 3301616701375337755ull},
    { "msf2.mtx", 8, 11227097946539390519ull},
    { "msf3.mtx", 5, 5965767602828141907ull},
    { "olm1000.mtx", 1, 14473458856538426155ull},
    { "pushpull.mtx", 1, 14764381483900318255ull},
    { "sample2.mtx", 1, 4918287057298807835ull},
    { "sample.mtx", 8, 11227097946539390519ull},
    { "structure.mtx", 3, 848279640410529436ull},
    { "test_BF.mtx", 3, 848279640410529436ull},
    { "test_FW_1000.mtx", 1, 14473458856538426155ull},
    { "test_FW_2003.mtx", 485, 13924889949050000093ull},
    { "test_FW_2500.mtx", 646, 9579946550331191330ull},
    { "tree-example.mtx", 1, 16959292359894689422ull},
    { "west0067_jumbled.mtx", 1, 15563611237648677666ull},
    { "west0067.mtx", 1, 15563611237648677666ull},
    { "west0067_noheader.mtx", 1, 15563611237648677666ull},
    { "zenios.mtx", 1391, 10773678236411609506ull},
    { "", 0, 0},
} ;
//------------------------------------------------------------------------------
// count_connected_components: count the # of components in a component vector
//------------------------------------------------------------------------------

int count_connected_components (GrB_Vector C, uint64_t *vector_hash) ;

int count_connected_components (GrB_Vector C, uint64_t *vector_hash)
{
    GrB_Index n = 0 ;
    OK (GrB_Vector_size (&n, C)) ;
    int ncomponents = 0 ;
    for (int i = 0 ; i < n ; i++)
    {
        int64_t comp = -1 ;
        int result = GrB_Vector_extractElement (&comp, C, i) ;
        if (result == GrB_SUCCESS && comp == i) ncomponents++ ;
        //hash all of the values into one number
        if (result == GrB_SUCCESS && vector_hash) 
        {
            (*vector_hash) *= 89734512321ull;
            (*vector_hash) += comp + i;
        }
    }
    return (ncomponents) ;
}

//****************************************************************************
void test_scc (void)
{
    #if LAGRAPH_SUITESPARSE
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
            uint64_t hash = 9238018047ull;
            int result_cc_count = count_connected_components(c, &hash);
            TEST_CHECK(result_cc_count == files[k].cc_count);
            TEST_CHECK(hash == files[k].hash);
            OK (LAGraph_Vector_Print (c, pr, stdout, msg)) ;
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
    #if LAGRAPH_SUITESPARSE
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
