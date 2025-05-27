//------------------------------------------------------------------------------
// LAGraph/experimental/test/test_diameter.c: test diameter
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

// TODO: vector results (peripheral, eccentricity, level, and parent) are not
// checked with exact values

#include <stdio.h>
#include <acutest.h>

#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Matrix C = NULL ;
GrB_Vector peripheral = NULL ;
GrB_Matrix level = NULL ;
GrB_Matrix parent = NULL ;
GrB_Vector sources = NULL ;
GrB_Vector est_peripheral = NULL ;
GrB_Vector eccentricity = NULL ;
GrB_Index n ;
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    int diameter ;
    bool symmetric ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    {  2, 1, "A.mtx" },
    { 60, 1, "jagmesh7.mtx" },
    {  6, 0, "west0067.mtx" }, // unsymmetric
    { 11, 1, "bcsstk13.mtx" },
    {  5, 1, "karate.mtx" },
    {  4, 1, "ldbc-cdlp-undirected-example.mtx" },
    {  4, 1, "ldbc-undirected-example-bool.mtx" },
    {  4, 1, "ldbc-undirected-example-unweighted.mtx" },
    {  4, 1, "ldbc-undirected-example.mtx" },
    {  3, 1, "ldbc-wcc-example.mtx" },
    {  0, 0, "" },
} ;

#undef OK
#define OK(method) \
{ \
    GrB_Info info = method ; \
    if (info != GrB_SUCCESS) \
    { \
        printf ("info: %d, msg: %s\n", info, msg) ; \
        TEST_CHECK (false) ; \
    } \
}

//------------------------------------------------------------------------------
// test_diameter
//------------------------------------------------------------------------------

void test_diameter (void)
{
    #if LAGRAPH_SUITESPARSE
    OK (LAGraph_Init (msg)) ;

    for (int k = 0 ; ; k++)
    {

        // load the matrix as A
        const char *aname = files [k].name ;
        bool symmetric = files [k].symmetric ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;
        OK (GrB_Matrix_nrows (&n, A)) ;
        LAGraph_PrintLevel pr = (n <= 100) ? LAGraph_COMPLETE : LAGraph_SHORT ;
        OK (LAGraph_Matrix_Print (A, pr, stdout, msg)) ;

        // construct a directed graph G with adjacency matrix S
        int kind = symmetric ?
            LAGraph_ADJACENCY_UNDIRECTED :
            LAGraph_ADJACENCY_DIRECTED ;
        OK (LAGraph_New (&G, &A, kind, msg)) ;
        TEST_CHECK (A == NULL) ;

        // compute the exact diameter
        GrB_Index diameter = 0 ;
        OK (LAGraph_ExactDiameter (&diameter, &peripheral, &eccentricity,
            G, 8, msg)) ;

        printf ("\n%s exact diameter: %" PRIu64 "\n", aname, diameter) ;
        TEST_CHECK (diameter == files [k].diameter) ;

        printf ("\nperipheral:\n") ;
        OK (LAGraph_Vector_Print (peripheral, pr, stdout, msg)) ;
        printf ("\neccentricity:\n") ;
        OK (LAGraph_Vector_Print (eccentricity, pr, stdout, msg)) ;

        for (int jit = 0 ; jit <= 1 ; jit++)
        {
            OK (LG_SET_JIT (jit ? GxB_JIT_ON : GxB_JIT_OFF)) ;
            // compute the estimated diameter
            GrB_Index estimated_diameter = 0 ;
            OK (LAGraph_EstimateDiameter (&estimated_diameter, &est_peripheral,
                G, 8, 4, 42, msg)) ;
            printf ("\nest diameter: %" PRIu64 "\n", estimated_diameter) ;
            TEST_CHECK (estimated_diameter <= diameter) ;
            printf ("\nest_peripheral:\n") ;
            OK (LAGraph_Vector_Print (est_peripheral, pr, stdout, msg)) ;
            OK (GrB_free (&est_peripheral)) ;
        }

        // try the multisource BFS directly
        OK (GrB_Vector_new (&sources, GrB_INT64, 4)) ;
        for (int i = 0 ; i < 4 ; i++)
        {
            OK (GrB_Vector_setElement (sources, i, i)) ;
        }
        OK (GrB_wait (sources, GrB_MATERIALIZE)) ;
        printf ("\nsource nodes for multiBFS:\n") ;
        OK (LAGraph_Vector_Print (sources, 5, stdout, msg)) ;
        OK (LAGraph_MultiSourceBFS (&level, &parent, G, sources, msg)) ;
        printf ("\nlevel:\n") ;
        OK (LAGraph_Matrix_Print (level, pr, stdout, msg)) ;
        printf ("\nparent:\n") ;
        OK (LAGraph_Matrix_Print (parent, pr, stdout, msg)) ;

        OK (GrB_free (&sources)) ;
        OK (GrB_free (&level)) ;
        OK (GrB_free (&parent)) ;
        OK (GrB_free (&peripheral)) ;
        OK (GrB_free (&eccentricity)) ;
        OK (LAGraph_Delete (&G, msg)) ;
    }

    OK (LAGraph_Finalize (msg)) ;
    #endif
}

//------------------------------------------------------------------------------
// test_diameter_huge
//------------------------------------------------------------------------------

void test_diameter_huge (void)
{
    #if LAGRAPH_SUITESPARSE
    OK (LAGraph_Init (msg)) ;
    OK (LG_SET_JIT (LG_JIT_OFF)) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&C, f, msg)) ;
    fclose (f) ;
    OK (GrB_Matrix_nrows (&n, C)) ;
    OK (LAGraph_Matrix_Print (C, 5, stdout, msg)) ;

    GrB_Index *I = NULL ;
    OK (LAGraph_Malloc ((void **) &I, n, sizeof (uint64_t), msg)) ;
    for (int k = 0 ; k < n ; k++)
    {
        I [k] = k ;
    }

    // A (0:n-1, 0:n-1) = C, where C is n-by-n
    OK (GrB_Matrix_new (&A, GrB_BOOL, UINT32_MAX, UINT32_MAX)) ;
    OK (GrB_assign (A, NULL, NULL, C, I, n, I, n, NULL)) ;

    // construct the graph
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    // compute the estimated diameter
    GrB_Index estimated_diameter = 0 ;
    OK (LAGraph_EstimateDiameter (&estimated_diameter, &est_peripheral,
        G, 8, 4, 42, msg)) ;
    printf ("\nest diameter: %" PRIu64 "\n", estimated_diameter) ;
    printf ("\nest_peripheral:\n") ;
    OK (LAGraph_Vector_Print (est_peripheral, 2, stdout, msg)) ;
    OK (GrB_free (&est_peripheral)) ;

    // try the multisource BFS directly
    OK (GrB_Vector_new (&sources, GrB_INT64, 4)) ;
    for (int i = 0 ; i < 4 ; i++)
    {
        OK (GrB_Vector_setElement (sources, i, i)) ;
    }
    OK (GrB_wait (sources, GrB_MATERIALIZE)) ;
    printf ("\nsource nodes for multiBFS:\n") ;
    OK (LAGraph_Vector_Print (sources, 5, stdout, msg)) ;
    OK (LAGraph_MultiSourceBFS (&level, &parent, G, sources, msg)) ;
    printf ("\nlevel:\n") ;
    OK (LAGraph_Matrix_Print (level, 2, stdout, msg)) ;
    printf ("\nparent:\n") ;
    OK (LAGraph_Matrix_Print (parent, 2, stdout, msg)) ;

    OK (GrB_free (&sources)) ;
    OK (GrB_free (&level)) ;
    OK (GrB_free (&parent)) ;
    OK (GrB_free (&C)) ;
    OK (LAGraph_Free ((void **) &I, msg)) ;
    OK (LAGraph_Delete (&G, msg)) ;

    OK (LAGraph_Finalize (msg)) ;
    #endif
}

//------------------------------------------------------------------------------
// test_errors
//------------------------------------------------------------------------------

void test_errors (void)
{
    #if LAGRAPH_SUITESPARSE
    LAGraph_Init (msg) ;
    // TODO: add error tests here
    LAGraph_Finalize (msg) ;
    #endif
}

//****************************************************************************

TEST_LIST = {
    {"diameter", test_diameter},
    {"diameter_huge", test_diameter_huge},
    {"diameter_errors", test_errors},
    {NULL, NULL}
} ;

