//------------------------------------------------------------------------------
// test_leiden.c: tests for LAGraph_Leiden community detection
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

#include <stdio.h>
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

typedef struct
{
    const char *matrix_file ;
    bool        expect_communities ;    // true iff non-trivial Q is expected
    double      min_modularity ;        // minimum acceptable modularity
    double      min_coverage ;          // minimum intra-cluster edge ratio
    double      min_performance ;       // minimum partition performance
} matrix_info ;

const matrix_info files[] =
{
    // matrix_file,      expect_communities, min_Q,  min_coverage, min_performance
    { "karate.mtx",      true,              0.35,   0.15,         0.30 },
    { "comm0.mtx",       true,              0.25,   0.10,         0.20 },
    { "",                false,             -1.0,   -1.0,         -1.0 }
} ;

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

        OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
        TEST_CHECK (A == NULL) ;    // LAGraph_New takes ownership

        // Ensure symmetry cache is populated (required by some checks).
        OK (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
        OK (LAGraph_Cached_OutDegree (G, msg)) ;

        uint64_t seed = 0 ; //unused
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
        int64_t min_label = 0, max_label = 0 ;

        OK (GrB_Vector_reduce_INT64 (
            &min_label, NULL, GxB_MIN_INT64_MONOID, c, NULL)) ;
        OK (GrB_Vector_reduce_INT64 (
            &max_label, NULL, GxB_MAX_INT64_MONOID, c, NULL)) ;

        TEST_CHECK (min_label >= 0) ;
        TEST_CHECK (max_label < n) ;

        GrB_Index n_communities = max_label + 1 ;

        // Compute modularity Q (requires SuiteSparse:GraphBLAS).
        double Q = 0.0 ;
        OK (LAGr_Modularity (&Q, 1.0, c, G, msg)) ;
        printf ("  Modularity Q = %f\n", Q) ;

        if (files[k].expect_communities)
        {
            // Validate modularity meets threshold for this graph
            TEST_CHECK (Q > files[k].min_modularity) ;
            TEST_MSG ("Expected Q > %f for %s, got Q = %f", 
                      files[k].min_modularity, aname, Q) ;
        }

        // Compute per-community statistics for robustness validation
        printf ("  Community statistics:\n") ;
        GrB_Index total_community_edges = 0 ;
        GrB_Index total_edges = 0 ;
        OK (GrB_Matrix_nvals (&total_edges, G->A)) ;

        for (GrB_Index com = 0 ; com < n_communities ; com++)
        {
            GrB_Index com_size = 0 ;
            GrB_Index com_edges = 0 ;

            // Count nodes and edges in this community
            for (GrB_Index i = 0 ; i < n ; i++)
            {
                int64_t label ;
                GrB_Info info = GrB_Vector_extractElement_INT64 (&label, c, i) ;
                if (info == GrB_SUCCESS && label == (int64_t)com)
                {
                    com_size++ ;
                    // Count internal edges: edges from i to other nodes in same community
                    for (GrB_Index j = 0 ; j < n ; j++)
                    {
                        int64_t label_j ;
                        GrB_Info info_j = GrB_Vector_extractElement_INT64 (&label_j, c, j) ;
                        if (info_j == GrB_SUCCESS && label_j == (int64_t)com)
                        {
                            bool has_edge = false ;
                            GrB_Info info_e = GrB_Matrix_extractElement_BOOL (&has_edge, G->A, i, j) ;
                            if (info_e == GrB_SUCCESS && has_edge)
                            {
                                com_edges++ ;
                            }
                        }
                    }
                }
            }

            if (com_size > 0)
            {
                double density = (com_size > 1) ? 
                    (double)com_edges / (com_size * (com_size - 1)) : 0.0 ;
                printf ("    Community %llu: size=%llu, edges=%llu, density=%f\n",
                        (unsigned long long) com, (unsigned long long) com_size,
                        (unsigned long long) com_edges, density) ;
                total_community_edges += com_edges ;
            }
        }

        // Compute edge-cut ratio (edges crossing communities / total edges)
        GrB_Index cut_edges = total_edges - total_community_edges ;
        double cut_ratio = (total_edges > 0) ? 
            (double)cut_edges / total_edges : 0.0 ;
        printf ("  Edge-cut ratio = %f (cut=%llu, total=%llu)\n", 
                cut_ratio, (unsigned long long) cut_edges, 
                (unsigned long long) total_edges) ;

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
