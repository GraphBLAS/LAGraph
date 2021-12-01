//------------------------------------------------------------------------------
// LAGraph/src/test/test_export.c: test LG_check_export
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include <LAGraph_test.h>

char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

#define LEN 512
char filename [LEN+1] ;

typedef struct
{
    LAGraph_Kind kind ;
    const char *name ;
}
matrix_info ;

const matrix_info files [ ] =
{
    { LAGRAPH_ADJACENCY_UNDIRECTED, "A.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "cover.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "jagmesh7.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "ldbc-cdlp-directed-example.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-cdlp-undirected-example.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "ldbc-directed-example.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-undirected-example.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "ldbc-wcc-example.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "LFAT5.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "msf1.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "msf2.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "msf3.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "sample2.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "sample.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "olm1000.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "bcsstk13.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "cryg2500.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "tree-example.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "west0067.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "karate.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_bool.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_int8.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_int16.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_int32.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_int64.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_uint8.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_uint16.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_int32.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "matrix_uint64.mtx" },
    { LAGRAPH_ADJACENCY_DIRECTED,   "skew_fp32.mtx" },
    { LAGRAPH_ADJACENCY_UNDIRECTED, "pushpull.mtx" },
    { LAGRAPH_UNKNOWN, "" },
} ;

//------------------------------------------------------------------------------
// test_export
//------------------------------------------------------------------------------

void test_export (void)
{
    LAGraph_Init (msg);
    GrB_Matrix A = NULL, C = NULL ;
    GrB_Type atype = NULL ;

    for (int k = 0 ; ; k++)
    {

        // load the adjacency matrix as A
        const char *aname = files [k].name ;
        LAGraph_Kind kind = files [k].kind ;
        if (strlen (aname) == 0) break;
        TEST_CASE (aname) ;
        printf ("\nMatrix: %s\n", aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;
        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        OK (LAGraph_MMRead (&A, &atype, f, msg)) ;
        OK (fclose (f)) ;
        TEST_MSG ("Loading of adjacency matrix failed") ;

        // create the graph
        OK (LAGraph_New (&G, &A, atype, kind, msg)) ;
        TEST_CHECK (A == NULL) ;    // A has been moved into G->A

        // export the graph
        GrB_Index *Ap = NULL ;
        GrB_Index *Aj = NULL ;
        void *Ax = NULL ;
        GrB_Index Ap_len, Aj_len, Ax_len ;
        size_t typesize ;

        OK (LG_check_export (G, &Ap, &Aj, &Ax, &Ap_len, &Aj_len,
            &Ax_len, &typesize, msg)) ;

        #if LG_SUITESPARSE
        OK (GxB_Matrix_import_CSR (&C
        #else
        LAGraph_Free ((void **) &Ap) ;
        LAGraph_Free ((void **) &Aj) ;
        LAGraph_Free ((void **) &Ax) ;
        #endif

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize(msg);
}

//****************************************************************************
//****************************************************************************
TEST_LIST = {
    {"test_export", test_export },
    {NULL, NULL}
};

