//------------------------------------------------------------------------------
// test_leiden.c: tests for LAGraph_Leiden community detection
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

#include <stdio.h>
#include <math.h>
#include <acutest.h>

#include "GraphBLAS.h"
#include "LG_Xtest.h"
#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg[LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
GrB_Matrix C = NULL ;
GrB_Scalar zero_bool = NULL ;
#define LEN 512
char filename[LEN + 1] ;

typedef struct
{
    const char *matrix_file ;
    LAGraph_Kind kind ;
    bool        force_symmetric ;       // make A symmetric via A + A'
    bool        force_positive ;        // make all edge weights positive
    bool        require_fp64 ;          // assert resulting matrix is FP64
    double      min_modularity ;        // set < -1 to skip threshold
} matrix_info ;

const matrix_info files[] =
{
    // matrix_file,      kind,                        force_sym, force_positive, require_fp64, min_Q
    { "karate.mtx",      LAGraph_ADJACENCY_UNDIRECTED, false,     true,          false,      0.35 },
    { "comm0.mtx",       LAGraph_ADJACENCY_UNDIRECTED, false,     true,          false,      0.25 },
    { "west0067.mtx",    LAGraph_ADJACENCY_DIRECTED,   true,      true,          true,      -2.0 },
    { "jagmesh7.mtx",    LAGraph_ADJACENCY_UNDIRECTED, false,     true,          false,     -2.0 },
    { "bcsstk13.mtx",    LAGraph_ADJACENCY_UNDIRECTED, false,     true,          true,      -2.0 },
    { "cryg2500.mtx",    LAGraph_ADJACENCY_DIRECTED,   true,      true,          true,      -2.0 },
    { "",                LAGraph_ADJACENCY_UNDIRECTED, false,     false,         false,     -2.0 }
} ;

const char *nonfinite_files[] =
{
    "matrix_fp64.mtx",
    "skew_fp64.mtx",
    ""
} ;

static int check_all_finite (const GrB_Matrix A, bool *all_finite)
{
    GxB_Iterator it = NULL ;
    GrB_Info info ;
    *all_finite = true ;
    info = GxB_Iterator_new (&it) ;
    if (info != GrB_SUCCESS) return info ;
    info = GxB_Matrix_Iterator_attach (it, A, NULL) ;
    if (info != GrB_SUCCESS) goto done ;
    info = GxB_Matrix_Iterator_seek (it, 0) ;
    while (info == GrB_SUCCESS)
    {
        double aij = GxB_Iterator_get_FP64 (it) ;
        if (!isfinite (aij))
        {
            *all_finite = false ;
            break ;
        }
        info = GxB_Matrix_Iterator_next (it) ;
    }
done:
    GrB_free (&it) ;
    if (info == GxB_EXHAUSTED || info == GrB_SUCCESS) return GrB_SUCCESS ;
    return info ;
}

//------------------------------------------------------------------------------
// test_Leiden
//------------------------------------------------------------------------------

void test_Leiden (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ;; k++)
    {
        const char *aname = files[k].matrix_file ;
        if (strlen (aname) == 0) break ;

        printf ("\n====== %s ======\n", aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;

        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        TEST_MSG ("Cannot open %s", filename) ;

        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;

        OK (LAGraph_New (&G, &A, files[k].kind, msg)) ;
        TEST_CHECK (A == NULL) ;    // LAGraph_New takes ownership

        // Ensure symmetry as needed for Leiden.
        if (files[k].force_symmetric)
        {
            OK (LAGraph_Cached_AT (G, msg)) ;
            OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
            if (G->is_symmetric_structure == LAGraph_FALSE)
            {
                OK (GrB_eWiseAdd (G->A, NULL, NULL, GrB_ONEB_FP64, G->A, G->AT, NULL)) ;
                G->is_symmetric_structure = LAGraph_TRUE ;
            }
        }
        else
        {
            OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
        }
        TEST_CHECK (G->is_symmetric_structure == LAGraph_TRUE) ;

        if (files[k].force_positive)
        {
            OK (GrB_apply (G->A, NULL, NULL, GrB_ABS_FP64, G->A, NULL)) ;
        }

        if (files[k].require_fp64)
        {
            GrB_Type atype = NULL ;
            OK (GxB_Matrix_type (&atype, G->A)) ;
            TEST_CHECK (atype == GrB_FP64) ;
        }

        OK (LAGraph_Cached_EMin (G, msg)) ;
        double min_val = 0.0 ;
        OK (GrB_Scalar_extractElement_FP64 (&min_val, G->emin)) ;
        TEST_CHECK (min_val >= 0.0) ;
        OK (LAGraph_Cached_OutDegree (G, msg)) ;

        uint64_t seed = 0 ; //unused
        GrB_Vector c = NULL ;

        GrB_Info info = LAGraph_Leiden (&c, G, seed, msg) ;
        TEST_CHECK (info == GrB_SUCCESS) ;
        if (info != GrB_SUCCESS)
        {
            GrB_free (&c) ;
            OK (LAGraph_Delete (&G, msg)) ;
            continue ;
        }
        TEST_CHECK (c != NULL) ;

        // Every node must have a community label.
        GrB_Index n, nvals ;
        OK (GrB_Matrix_nrows (&n, G->A)) ;
        OK (GrB_Vector_nvals (&nvals, c)) ;
        TEST_CHECK (nvals == n) ;
        TEST_MSG ("Expected all %llu nodes to have labels, got %llu",
                  (unsigned long long) n, (unsigned long long) nvals) ;

        // Community labels must be in [0, n-1].
        int64_t min_label = 0, max_label = 0 ;

        OK (GrB_Vector_reduce_INT64 (
            &min_label, NULL, GxB_MIN_INT64_MONOID, c, NULL)) ;
        OK (GrB_Vector_reduce_INT64 (
            &max_label, NULL, GxB_MAX_INT64_MONOID, c, NULL)) ;

        TEST_CHECK (min_label >= 0) ;
        TEST_CHECK (max_label < n) ;

        // Compute modularity Q
        double Q = 0.0 ;
        GrB_Descriptor desc = NULL;
        OK (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;
        OK (GrB_Scalar_new (&zero_bool, GrB_BOOL)) ;
        OK (GrB_Scalar_setElement_BOOL (zero_bool, false)) ;
        OK (GrB_Descriptor_new (&desc)) ;
        OK (GrB_set (desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)) ;
        OK (GxB_Matrix_build_Scalar_Vector(C, c, c, zero_bool, desc)) ;
        OK (LAGr_AdjModularity (&Q, 1.0, G->A, C, msg)) ;
        printf ("  Modularity Q = %f\n", Q) ;

        if (files[k].min_modularity > -1.0)
        {
            // Validate modularity meets threshold for this graph
            TEST_CHECK (Q > files[k].min_modularity) ;
            TEST_MSG ("Expected Q > %f for %s, got Q = %f", 
                      files[k].min_modularity, aname, Q) ;
        }
        TEST_CHECK (isfinite (Q)) ;
        TEST_CHECK (Q >= -1.0 && Q <= 1.0) ;

        GrB_free (&c) ;
        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize (msg) ;
}

void test_Leiden_NonfiniteInputs (void)
{
    LAGraph_Init (msg) ;

    for (int k = 0 ;; k++)
    {
        const char *aname = nonfinite_files [k] ;
        if (strlen (aname) == 0) break ;

        printf ("\n====== nonfinite %s ======\n", aname) ;
        snprintf (filename, LEN, LG_DATA_DIR "%s", aname) ;

        FILE *f = fopen (filename, "r") ;
        TEST_CHECK (f != NULL) ;
        TEST_MSG ("Cannot open %s", filename) ;
        OK (LAGraph_MMRead (&A, f, msg)) ;
        fclose (f) ;

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;
        OK (LAGraph_Cached_AT (G, msg)) ;
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
        if (G->is_symmetric_structure == LAGraph_FALSE)
        {
            OK (GrB_eWiseAdd (G->A, NULL, NULL, GrB_PLUS_FP64, G->A, G->AT, NULL)) ;
            G->is_symmetric_structure = LAGraph_TRUE ;
        }

        GrB_Type atype = NULL ;
        OK (GxB_Matrix_type (&atype, G->A)) ;
        TEST_CHECK (atype == GrB_FP64) ;

        bool all_finite = true ;
        OK (check_all_finite (G->A, &all_finite)) ;
        TEST_CHECK (!all_finite) ;

        OK (LAGraph_Delete (&G, msg)) ;
    }

    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test list
//------------------------------------------------------------------------------

TEST_LIST =
{
    { "Leiden", test_Leiden },
    { "Leiden nonfinite inputs", test_Leiden_NonfiniteInputs },
    { NULL, NULL }
} ;
