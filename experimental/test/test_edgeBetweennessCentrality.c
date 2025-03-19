//------------------------------------------------------------------------------
// LAGraph/src/test/test_Betweenness.c: test cases for BC (GAP method)
// -----------------------------------------------------------------------------

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

#include <stdio.h>
#include <acutest.h>
#include "LAGraphX.h" // important
#include "LAGraph_test.h"
#include "LG_Xtest.h"
#include "LG_internal.h"

#define LEN 512
char msg [LAGRAPH_MSG_LEN] ;
char filename [LEN+1] ;
LAGraph_Graph G = NULL ;

//------------------------------------------------------------------------------
// difference: compare the LAGraph and GAP results
//------------------------------------------------------------------------------

double difference(GrB_Matrix bc, double* reference_bc, GrB_Index rows, GrB_Index cols) ;

double difference(GrB_Matrix bc, double* reference_bc, GrB_Index rows, GrB_Index cols)
{
    // GrB_Matrix diff = NULL;
    GrB_Matrix diff = NULL, reference_bc_matrix = NULL ;
    OK(GrB_Matrix_new(&reference_bc_matrix, GrB_FP64, rows, cols)) ;

    // Populate gap_bc with values from gap_result
    for (GrB_Index i = 0; i < rows; i++) {
        for (GrB_Index j = 0; j < cols; j++) {
            OK(GrB_Matrix_setElement_FP64(reference_bc_matrix, *(reference_bc + i * cols + j), i, j)) ;
        }
    }

    // Compute diff = max(abs(reference_bc_matrix - bc))
    OK(GrB_Matrix_new(&diff, GrB_FP64, rows, cols)) ;
    OK(GrB_eWiseAdd(diff, NULL, NULL, GrB_MINUS_FP64, reference_bc_matrix, bc, NULL)) ;
    OK(GrB_apply(diff, NULL, NULL, GrB_ABS_FP64, diff, NULL)) ;

    double err = 0 ;
    OK(GrB_reduce(&err, NULL, GrB_MAX_MONOID_FP64, diff, NULL)) ;

    OK(GrB_free(&diff)) ;
    OK(GrB_free(&reference_bc_matrix)) ;

    return err ;
}

double matrix_difference(GrB_Matrix bc, GrB_Matrix reference_bc) ;

double matrix_difference(GrB_Matrix bc, GrB_Matrix reference_bc)
{
    GrB_Matrix diff = NULL ;

    uint64_t n ;
    GrB_Matrix_nrows (&n, bc) ;

    // Compute diff = max(abs(reference_bc - bc))
    GrB_Matrix_new(&diff, GrB_FP64, n, n) ;
    GrB_eWiseAdd(diff, NULL, NULL, GrB_MINUS_FP64, reference_bc, bc, NULL) ;
    GrB_apply(diff, NULL, NULL, GrB_ABS_FP64, diff, NULL) ;

    double err = 1 ;
    GrB_reduce(&err, NULL, GrB_MAX_MONOID_FP64, diff, NULL) ;

    GrB_free(&diff) ;

    return err ;
}

//------------------------------------------------------------------------------
// results for book graph
//------------------------------------------------------------------------------

// Exact results from edge_betweenness_centrality from NetworkX of the graph
// from the book

double diamonds_ebc [8][8] = 
{
    {0.0, 2.333333333333333, 2.333333333333333, 2.333333333333333, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0, 5.333333333333333, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 5.333333333333333, 0.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0, 5.333333333333333, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 7.5, 7.5, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.5},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
} ; 

//------------------------------------------------------------------------------
// results for karate graph
//------------------------------------------------------------------------------

// Exact results from edge_betweenness_centrality from NetworkX of the karate 
// graph

double karate_ebc [34][34] = 
{
    {0.0, 14.166666666666664, 43.638888888888886, 11.5, 29.333333333333332, 43.83333333333333, 43.833333333333336, 12.80238095238095, 41.64841269841271, 0.0, 29.333333333333332, 33.0, 26.099999999999994, 23.77063492063493, 0.0, 0.0, 0.0, 22.509523809523813, 0.0, 25.770634920634926, 0.0, 22.50952380952381, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 71.39285714285712, 0.0, 0.0},
    {14.166666666666664, 0.0, 13.033333333333335, 4.333333333333333, 0.0, 0.0, 0.0, 4.164285714285714, 0.0, 0.0, 0.0, 0.0, 0.0, 6.9595238095238106, 0.0, 0.0, 0.0, 10.490476190476187, 0.0, 8.209523809523809, 0.0, 10.490476190476187, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.10952380952381, 0.0, 0.0, 0.0},
    {43.638888888888886, 13.033333333333335, 0.0, 12.583333333333332, 0.0, 0.0, 0.0, 14.145238095238092, 5.147619047619047, 17.28095238095238, 0.0, 0.0, 0.0, 4.28095238095238, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.10873015873016, 12.780952380952376, 0.0, 0.0, 0.0, 38.70158730158729, 0.0},
    {11.5, 4.333333333333333, 12.583333333333332, 0.0, 0.0, 0.0, 0.0, 1.8880952380952383, 0.0, 0.0, 0.0, 0.0, 6.899999999999997, 8.37142857142857, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {29.333333333333332, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6666666666666665, 0.0, 0.0, 0.0, 1.6666666666666665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {43.83333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6666666666666667, 0.0, 0.0, 0.0, 2.6666666666666665, 0.0, 0.0, 0.0, 0.0, 0.0, 16.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {43.833333333333336, 0.0, 0.0, 0.0, 2.6666666666666665, 1.6666666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {12.80238095238095, 4.164285714285714, 14.145238095238092, 1.8880952380952383, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {41.64841269841271, 0.0, 5.147619047619047, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.5, 0.0, 17.077777777777776, 22.684920634920633},
    {0.0, 0.0, 17.28095238095238, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.614285714285714},
    {29.333333333333332, 0.0, 0.0, 0.0, 1.6666666666666665, 2.6666666666666665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {26.099999999999994, 0.0, 0.0, 6.899999999999997, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {23.77063492063493, 6.9595238095238106, 4.28095238095238, 8.37142857142857, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 38.04920634920634},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.511111111111113, 19.488888888888887},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.511111111111113, 19.488888888888887},
    {0.0, 0.0, 0.0, 0.0, 0.0, 16.5, 16.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {22.509523809523813, 10.490476190476187, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.511111111111113, 19.488888888888887},
    {25.770634920634926, 8.209523809523809, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 33.31349206349207},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.511111111111113, 19.488888888888887},
    {22.50952380952381, 10.490476190476187, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.511111111111111, 19.488888888888887},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.094444444444443, 0.0, 5.9111111111111105, 0.0, 3.7333333333333334, 0.0, 0.0, 12.533333333333331, 18.327777777777783},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3666666666666667, 0.0, 10.466666666666665, 0.0, 0.0, 0.0, 22.5, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.094444444444443, 2.3666666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.59444444444445, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.5428571428571427, 0.0, 0.0, 0.0, 30.457142857142856},
    {0.0, 0.0, 23.10873015873016, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.9111111111111105, 10.466666666666665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.097619047619048},
    {0.0, 0.0, 12.780952380952376, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.333333333333332, 0.0, 13.78095238095238},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.7333333333333334, 0.0, 0.0, 2.5428571428571427, 0.0, 0.0, 0.0, 0.0, 0.0, 13.087301587301585, 16.72222222222222},
    {0.0, 18.10952380952381, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.566666666666666, 15.042857142857141},
    {71.39285714285712, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.5, 23.59444444444445, 0.0, 0.0, 8.333333333333332, 0.0, 0.0, 0.0, 23.244444444444447, 29.95396825396826},
    {0.0, 0.0, 38.70158730158729, 0.0, 0.0, 0.0, 0.0, 0.0, 17.077777777777776, 0.0, 0.0, 0.0, 0.0, 0.0, 13.511111111111113, 13.511111111111113, 0.0, 0.0, 13.511111111111113, 0.0, 13.511111111111113, 0.0, 13.511111111111111, 12.533333333333331, 0.0, 0.0, 0.0, 0.0, 0.0, 13.087301587301585, 9.566666666666666, 23.244444444444447, 0.0, 4.614285714285714},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.684920634920633, 16.614285714285714, 0.0, 0.0, 0.0, 38.04920634920634, 19.488888888888887, 19.488888888888887, 0.0, 0.0, 19.488888888888887, 33.31349206349207, 19.488888888888887, 0.0, 19.488888888888887, 18.327777777777783, 0.0, 0.0, 30.457142857142856, 17.097619047619048, 13.78095238095238, 16.72222222222222, 15.042857142857141, 29.95396825396826, 4.614285714285714, 0.0},
} ; 


//------------------------------------------------------------------------------
// test_diamonds_ebc
//------------------------------------------------------------------------------

void test_diamonds_ebc (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Matrix AT = NULL ;
    GrB_Matrix centrality = NULL ;
    int niters = 0 ;
    LAGraph_Kind kind = LAGraph_ADJACENCY_DIRECTED ;

    // create the karate graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "diamonds.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, kind, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A

    // check that AT is cached
    int ok_result = (kind == LAGraph_ADJACENCY_UNDIRECTED) ?
        LAGRAPH_CACHE_NOT_NEEDED : GrB_SUCCESS ;
    int result = LAGraph_Cached_AT (G, msg) ;
    TEST_CHECK (result == ok_result) ;

    // compute its betweenness centrality with C version
    OK (LG_check_edgeBetweennessCentrality (&centrality, G, msg)) ;
    double err = difference(centrality, &diamonds_ebc[0][0], 8, 8) ;
    printf ("\n  diamonds:   err: %e (C version)", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    // compute its betweenness centrality with GraphBLAS version
    OK (LAGr_EdgeBetweennessCentrality (&centrality, G, msg)) ;
    err = difference(centrality, &diamonds_ebc[0][0], 8, 8) ;
    printf ("\n  diamonds:   err: %e (pure GraphBLAS)\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    OK (LAGraph_Delete (&G, msg)) ;
    LAGraph_Finalize (msg) ;
}

//------------------------------------------------------------------------------
// test_karate_ebc
//------------------------------------------------------------------------------

void test_karate_ebc (void)
{
    LAGraph_Init (msg) ;
    GrB_Matrix A = NULL ;
    GrB_Matrix centrality = NULL ;
    int niters = 0 ;

    // create the karate graph
    snprintf (filename, LEN, LG_DATA_DIR "%s", "karate.mtx") ;
    FILE *f = fopen (filename, "r") ;
    TEST_CHECK (f != NULL) ;
    OK (LAGraph_MMRead (&A, f, msg)) ;
    OK (fclose (f)) ;
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A

    // compute its betweenness centrality (C version)
    OK (LG_check_edgeBetweennessCentrality (&centrality, G, msg)) ;
    double err = difference(centrality, &karate_ebc[0][0], 34, 34) ;
    printf ("\n  karate:   err: %e (C version)", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    // compute its betweenness centrality (GraphBLAS version)
    OK (LAGr_EdgeBetweennessCentrality (&centrality, G, msg)) ;
    err = difference(centrality, &karate_ebc[0][0], 34, 34) ;
    printf ("\n  karate:   err: %e (GraphBLAS version)\n", err) ;
    TEST_CHECK (err < 1e-4) ;
    OK (GrB_free (&centrality)) ;

    OK (LAGraph_Delete (&G, msg)) ;
    LAGraph_Finalize (msg) ;

}

// Function to test multiple matrix market files
void test_many(void)
{
    LAGraph_Init(msg);

    const char *files[] = {
        "random_unweighted_bipartite1.mtx",
        "random_unweighted_bipartite2.mtx",
        "random_unweighted_general1.mtx",
        "random_unweighted_general2.mtx",
        NULL
    };

    for (int i = 0; files[i] != NULL; i++)
    {
        GrB_Matrix A = NULL;
        GrB_Matrix centrality = NULL;
        GrB_Matrix reference_centrality = NULL;

        snprintf(filename, LEN, LG_DATA_DIR "%s", files[i]);
        FILE *f = fopen(filename, "r");
        TEST_CHECK(f != NULL);
        OK(LAGraph_MMRead(&A, f, msg));
        OK(fclose(f));
        OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_DIRECTED, msg));
        LAGRAPH_TRY (LAGraph_DeleteSelfEdges (G, msg)) ;
        LAGRAPH_TRY (LAGraph_Cached_AT (G, msg)) ;
        TEST_CHECK(A == NULL); // A has been moved into G->A

        // compute its betweenness centrality (GraphBLAS version)
        OK(LAGr_EdgeBetweennessCentrality(&centrality, G, msg));

        // compute its betweenness centrality (C version)
        OK(LG_check_edgeBetweennessCentrality(&reference_centrality, G, msg));

        // Compare the results
        double err = matrix_difference(centrality, reference_centrality);
        printf("\n  %s: err: %e", files[i], err);
        TEST_CHECK(err < 1e-4);

        OK(GrB_free(&centrality));
        OK(GrB_free(&reference_centrality));
        OK(LAGraph_Delete(&G, msg));
    }
    printf("\n") ;

    LAGraph_Finalize(msg);
}

//------------------------------------------------------------------------------
// list of tests
//------------------------------------------------------------------------------


TEST_LIST = {
    {"test_diamonds_ebc", test_diamonds_ebc},
    {"test_karate_ebc", test_karate_ebc},
    {"test_many", test_many},
    {NULL, NULL}
};
