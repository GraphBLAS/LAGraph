//----------------------------------------------------------------------------
// LAGraph/experimental/test/test_AllKtest.c: test cases for all-k-truss
// ----------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//-----------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>

#include "LAGraphX.h"
#include "LAGraph_test.h"
#include "LG_Xtest.h"

char msg [LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Matrix C1 = NULL ;
GrB_Type atype = NULL ;
#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    uint32_t ntriangles ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    {     11, "A.mtx" },
    {   2016, "jagmesh7.mtx" },
    { 342300, "bcsstk13.mtx" },
    {     45, "karate.mtx" },
    {      6, "ldbc-cdlp-undirected-example.mtx" },
    {      4, "ldbc-undirected-example-bool.mtx" },
    {      4, "ldbc-undirected-example-unweighted.mtx" },
    {      4, "ldbc-undirected-example.mtx" },
    {      5, "ldbc-wcc-example.mtx" },
    { 0, "" },
} ;

//****************************************************************************
void test_AllKTruss (void)
{
    LAGraph_Init (msg) ;

    for (int id = 0 ; ; id++)
    {

        // load the matrix as A
        const char *aname = files [id].name ;
        uint32_t ntriangles = files [id].ntriangles ;
        if (strlen (aname) == 0) break;
        printf ("\n================================== %s:\n", aname) ;
        TEST_CASE (aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;
        fclose (f) ;

        // construct an undirected graph G with adjacency matrix A
        OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;

        // check for self-edges
        OK (LAGraph_Property_NDiag (G, msg)) ;
        if (G->ndiag != 0)
        {
            // remove self-edges
            printf ("graph has %g self edges\n", (double) G->ndiag) ;
            OK (LAGraph_DeleteDiag (G, msg)) ;
            printf ("now has %g self edges\n", (double) G->ndiag) ;
            TEST_CHECK (G->ndiag == 0) ;
        }

        // compute each k-truss
        bool ok = false ;
        GrB_Index n ;
        int64_t kmax ;
        OK (GrB_Matrix_nrows (&n, G->A)) ;
        GrB_Matrix *Cset = (GrB_Matrix *) LAGraph_Calloc (n,
            sizeof (GrB_Matrix)) ;
        int64_t *ntris  = LAGraph_Malloc (n, sizeof (int64_t)) ;
        int64_t *nedges = LAGraph_Malloc (n, sizeof (int64_t)) ;
        int64_t *nsteps = LAGraph_Malloc (n, sizeof (int64_t)) ;
        OK (LAGraph_AllKTruss (Cset, &kmax, ntris, nedges, nsteps, G, msg)) ;
        printf ("all k-truss: kmax %g\n", (double) kmax) ;

        // compute each k-truss using LAGraph_KTruss, and compare
        for (int k = 3 ; k < n ; k++)
        {
            // printf ("\n%d-truss:\n", k) ;
            TEST_CHECK (k <= kmax) ;
            // compute the k-truss
            OK (LAGraph_KTruss (&C1, G, k, msg)) ;

            // check the result
            GrB_Index nvals ;
            OK (GrB_Matrix_nvals (&nvals, C1)) ;
            OK (LAGraph_IsEqual (&ok, C1, Cset [k], msg)) ;
            TEST_CHECK (ok) ;

            // count the triangles in the 3-truss
            uint32_t nt = 0 ;
            OK (GrB_reduce (&nt, NULL, GrB_PLUS_MONOID_UINT32, C1, NULL)) ;
            nt = nt / 6 ;
            if (k == 3)
            {
                TEST_CHECK (nt == ntriangles) ;
            }
            TEST_CHECK (nt == ntris [k]) ;
            TEST_CHECK (nvals == 2 * nedges [k]) ;
            TEST_CHECK (nsteps [k] >= 0) ;

            // free C1, and break if C1 is empty
            OK (GrB_free (&C1)) ;
            if (nvals == 0)
            {
                TEST_CHECK (k == kmax) ;
                break ;
            }
        }

        // convert to directed with symmetric structure and recompute
        G->kind = LAGRAPH_ADJACENCY_DIRECTED ;
        G->A_structure_is_symmetric = true ;
        int64_t k2 ;
        GrB_Matrix *Cset2 = (GrB_Matrix *) LAGraph_Calloc (n,
            sizeof (GrB_Matrix)) ;
        int64_t *ntris2  = LAGraph_Malloc (n, sizeof (int64_t)) ;
        int64_t *nedges2 = LAGraph_Malloc (n, sizeof (int64_t)) ;
        int64_t *nsteps2 = LAGraph_Malloc (n, sizeof (int64_t)) ;
        OK (LAGraph_AllKTruss (Cset2, &k2, ntris2, nedges2, nsteps2, G, msg)) ;
        TEST_CHECK (k2 == kmax) ;
        for (int k = 0 ; k <= kmax ; k++)
        {
            TEST_CHECK (ntris2  [k] == ntris  [k]) ;
            TEST_CHECK (nedges2 [k] == nedges [k]) ;
            TEST_CHECK (nsteps2 [k] == nsteps [k]) ;
            if (k < 3)
            {
                TEST_CHECK (Cset [k] == NULL) ;
                TEST_CHECK (Cset2 [k] == NULL) ;
            }
            else
            {
                OK (LAGraph_IsEqual (&ok, Cset [k], Cset2 [k], msg)) ;
            }
//          if (!ok)
//          {
//              GxB_print (Cset [k], 2) ;
//              GxB_print (Cset2 [k], 2) ;
//          }
            TEST_CHECK (ok) ;
            OK (GrB_free (&(Cset [k]))) ;
            OK (GrB_free (&(Cset2 [k]))) ;
        }

        LAGraph_Free ((void **) &Cset) ;
        LAGraph_Free ((void **) &ntris) ;
        LAGraph_Free ((void **) &nedges) ;
        LAGraph_Free ((void **) &nsteps) ;

        LAGraph_Free ((void **) &Cset2) ;
        LAGraph_Free ((void **) &ntris2) ;
        LAGraph_Free ((void **) &nedges2) ;
        LAGraph_Free ((void **) &nsteps2) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_AllKTruss_errors
//------------------------------------------------------------------------------

void test_allktruss_errors (void)
{
    LAGraph_Init (msg) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
    TEST_MSG ("Loading of adjacency matrix failed") ;
    fclose (f) ;

    // construct an undirected graph G with adjacency matrix A
    OK (LAGraph_New (&G, &A, atype, LAGRAPH_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;

    OK (LAGraph_Property_NDiag (G, msg)) ;

    GrB_Index n ;
    int64_t kmax ;
    OK (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_Matrix *Cset = (GrB_Matrix *) LAGraph_Calloc (n, sizeof (GrB_Matrix)) ;
    int64_t *ntris  = LAGraph_Malloc (n, sizeof (int64_t)) ;
    int64_t *nedges = LAGraph_Malloc (n, sizeof (int64_t)) ;
    int64_t *nsteps = LAGraph_Malloc (n, sizeof (int64_t)) ;

    // kmax is NULL
    int result = LAGraph_AllKTruss (Cset, NULL, ntris, nedges, nsteps, G, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_NULL_POINTER) ;

    // G is invalid
    result = LAGraph_AllKTruss (Cset, &kmax, ntris, nedges, nsteps, NULL, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == GrB_INVALID_OBJECT) ;

    // G may have self edges
    G->ndiag = LAGRAPH_UNKNOWN ;
    result = LAGraph_AllKTruss (Cset, &kmax, ntris, nedges, nsteps, G, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == -1004) ;

    // G is undirected
    G->ndiag = 0 ;
    G->kind = LAGRAPH_ADJACENCY_DIRECTED ;
    G->A_structure_is_symmetric = LAGRAPH_FALSE ;
    result = LAGraph_AllKTruss (Cset, &kmax, ntris, nedges, nsteps, G, msg) ;
    printf ("\nresult: %d %s\n", result, msg) ;
    TEST_CHECK (result == -1005) ;

    LAGraph_Free ((void **) &Cset) ;
    LAGraph_Free ((void **) &ntris) ;
    LAGraph_Free ((void **) &nedges) ;
    LAGraph_Free ((void **) &nsteps) ;

    OK (LAGraph_Delete (&G, msg)) ;
    LAGraph_Finalize (msg) ;
}

//****************************************************************************

TEST_LIST = {
    {"allktruss", test_AllKTruss},
    {"allktruss_errors", test_allktruss_errors},
    {NULL, NULL}
};
