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
#define LEN 512
char filename[LEN + 1] ;

int LG_Leiden_move_nodes
(
    GrB_Matrix C,
    uint64_t *community,
    uint64_t *nodes_popped_handle,
    uint64_t *nodes_evaluated_handle,
    uint64_t *nodes_moved_handle,
    const GrB_Matrix A,
    const GrB_Vector deg,
    uint64_t *queue,
    bool *enqueued,
    double m_inv2,
    char *msg
) ;

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
        GrB_Matrix C = NULL ;
        GrB_Scalar zero_bool = NULL ;
        GrB_Descriptor desc = NULL ;
        OK (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;
        OK (GrB_Scalar_new (&zero_bool, GrB_BOOL)) ;
        OK (GrB_Scalar_setElement_BOOL (zero_bool, false)) ;
        OK (GrB_Descriptor_new (&desc)) ;
        OK (GrB_set (desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)) ;
        OK (GxB_Matrix_build_Scalar_Vector (C, c, c, zero_bool, desc)) ;
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

        GrB_free (&desc) ;
        GrB_free (&zero_bool) ;
        GrB_free (&C) ;
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

void test_Leiden_Phase1SingletonSplit (void)
{
    LAGraph_Init (msg) ;

    const GrB_Index n = 4 ;
    GrB_Matrix B = NULL, C = NULL ;
    GrB_Vector deg = NULL ;
    OK (GrB_Matrix_new (&B, GrB_FP64, n, n)) ;
    OK (GrB_Matrix_new (&C, GrB_BOOL, n, n)) ;
    OK (GrB_Vector_new (&deg, GrB_FP64, n)) ;

    GrB_Index Ai [4] = { 0, 2, 0, 3 } ;
    GrB_Index Aj [4] = { 2, 0, 3, 0 } ;
    double Ax [4] = { 1.0, 1.0, 1.0, 1.0 } ;
    OK (GrB_Matrix_build_FP64 (B, Ai, Aj, Ax, 4, GrB_PLUS_FP64)) ;

    uint64_t community [4] = { 0, 0, 2, 3 } ;
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        OK (GrB_Matrix_setElement_BOOL (C, true, i, community [i])) ;
    }

    // Deliberately chosen to force negative neighbor gains for node 0.
    OK (GrB_Vector_setElement_FP64 (deg, 100.0, 0)) ;
    OK (GrB_Vector_setElement_FP64 (deg, 100.0, 1)) ;
    OK (GrB_Vector_setElement_FP64 (deg, 100.0, 2)) ;
    OK (GrB_Vector_setElement_FP64 (deg, 100.0, 3)) ;

    uint64_t queue [5] = { 0, 1, 2, 3, 0 } ;
    bool enqueued [4] = { true, true, true, true } ;
    uint64_t popped = 0, evaluated = 0, moved = 0 ;
    OK (LG_Leiden_move_nodes (C, community, &popped, &evaluated, &moved,
        B, deg, queue, enqueued, -1.0 / 400.0, msg)) ;

    // Community 1 starts empty; node 0 should split into it.
    TEST_CHECK (community [0] == 1) ;

    OK (GrB_free (&deg)) ;
    OK (GrB_free (&C)) ;
    OK (GrB_free (&B)) ;
    LAGraph_Finalize (msg) ;
}

#if LG_BRUTAL_TESTS
void test_Leiden_brutal (void)
{
    OK (LG_brutal_setup (msg)) ;
    OK (GxB_Global_Option_set (GxB_JIT_C_CONTROL, GxB_JIT_OFF)) ;

    snprintf (filename, LEN, LG_DATA_DIR "%s", "comm0.mtx") ;
    uint64_t seed = 0 ;
    GrB_Vector c = NULL ;
    LAGraph_Graph H = NULL ;
    GrB_Matrix B = NULL ;

    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    TEST_MSG ("Cannot open %s", filename) ;
    OK (LAGraph_MMRead (&B, f, msg)) ;
    fclose (f) ;
    OK (LAGraph_New (&H, &B, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    OK (LAGraph_Cached_IsSymmetricStructure (H, msg)) ;
    TEST_CHECK (H->is_symmetric_structure == LAGraph_TRUE) ;
    OK (GrB_apply (H->A, NULL, NULL, GrB_ABS_FP64, H->A, NULL)) ;
    OK (LAGraph_Cached_EMin (H, msg)) ;

    LG_brutal = INT64_MAX ;
    OK (LAGraph_Leiden (&c, H, seed, msg)) ;
    LG_brutal = -1 ;

    GrB_Index n, nvals ;
    OK (GrB_Matrix_nrows (&n, H->A)) ;
    OK (GrB_Vector_nvals (&nvals, c)) ;
    TEST_CHECK (nvals == n) ;

    GrB_free (&c) ;
    OK (LAGraph_Delete (&H, msg)) ;
    GrB_free (&B) ;

    OK (LG_brutal_teardown (msg)) ;
}
#endif

//------------------------------------------------------------------------------
// test list
//------------------------------------------------------------------------------

TEST_LIST =
{
    { "Leiden", test_Leiden },
    { "Leiden nonfinite inputs", test_Leiden_NonfiniteInputs },
    { "Leiden phase1 singleton split", test_Leiden_Phase1SingletonSplit },
    #if LG_BRUTAL_TESTS
    { "Leiden brutal", test_Leiden_brutal },
    #endif
    { NULL, NULL }
} ;
