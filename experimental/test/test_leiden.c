//------------------------------------------------------------------------------
// test_leiden.c: tests for LAGraph_Leiden community detection
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

#include <stdio.h>
#include <acutest.h>

#include "LG_Xtest.h"
#include <LAGraphX.h>
#include <LAGraph_test.h>

char msg[LAGRAPH_MSG_LEN] ;
LAGraph_Graph G = NULL ;
GrB_Matrix A = NULL ;
#define LEN 512
char filename[LEN + 1] ;

typedef struct
{
    const char *matrix_file ;
    bool        expect_communities ;    // true iff non-trivial Q is expected
} matrix_info ;

const matrix_info files[] =
{
    { "karate.mtx", true  },    // Zachary karate club (Q > 0 expected)
    { "",           false }
} ;

//------------------------------------------------------------------------------
// test_Leiden
//------------------------------------------------------------------------------

void test_Leiden (void)
{
    LAGraph_Init (msg) ;

#if LAGRAPH_SUITESPARSE
    // Disable JIT before any GraphBLAS operations; the JIT compiler may not
    // be available in all build environments.
    OK (LG_SET_JIT (LG_JIT_OFF)) ;
#endif

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

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;    // LAGraph_New takes ownership

        // Ensure symmetry cache is populated (required by some checks).
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;

        uint64_t seed = 42 ;
        GrB_Vector c = NULL ;

        OK (LAGraph_Leiden (&c, G, seed, msg)) ;
        TEST_CHECK (c != NULL) ;

        // Every node must have a community label.
        GrB_Index n, nvals ;
        OK (GrB_Matrix_nrows (&n, G->A)) ;
        OK (GrB_Vector_nvals (&nvals, c)) ;
        TEST_CHECK (nvals == n) ;
        TEST_MSG ("Expected all %llu nodes to have labels, got %llu",
                  (unsigned long long) n, (unsigned long long) nvals) ;

        // Community labels must be in [0, n-1].
        for (GrB_Index i = 0 ; i < n ; i++)
        {
            int64_t label = 0 ;
            OK (GrB_Vector_extractElement_INT64 (&label, c, i)) ;
            TEST_CHECK (label >= 0 && (GrB_Index) label < n) ;
        }

#if LAGRAPH_SUITESPARSE
        // Compute modularity Q (requires SuiteSparse:GraphBLAS).
        double Q = 0.0 ;
        OK (LAGr_Modularity (&Q, 1.0, c, G, msg)) ;
        printf ("  Modularity Q = %f  (K_ref <= %llu)\n",
                Q, (unsigned long long) n) ;

        if (files[k].expect_communities)
        {
            // Multi-level Leiden on karate.mtx achieves Q ≈ 0.42, matching
            // published Louvain/Leiden benchmarks on this graph.
            TEST_CHECK (Q > 0.37) ;
            TEST_MSG ("Expected Q > 0.37 for %s, got Q = %f", aname, Q) ;
        }
#endif

        GrB_free (&c) ;
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
    { NULL, NULL }
} ;
