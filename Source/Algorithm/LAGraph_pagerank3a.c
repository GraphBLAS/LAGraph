//------------------------------------------------------------------------------
// LAGraph_pagerank3a: pagerank using a real semiring
//------------------------------------------------------------------------------

/*
LAGraph:  graph algorithms based on GraphBLAS

Copyright 2019 LAGraph Contributors.

(see Contributors.txt for a full list of Contributors; see
ContributionInstructions.txt for information on how you can Contribute to
this project).

All Rights Reserved.

NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

Released under a BSD license, please see the LICENSE file distributed with
this Software or contact permission@sei.cmu.edu for full terms.

Created, in part, with funding and support from the United States
Government.  (see Acknowledgments.txt file).

This program includes and/or can make use of certain third party source
code, object code, documentation and other files ("Third Party Software").
See LICENSE file for more details.

*/

// LAGraph_pagerank3a: Alternative PageRank implementation using a real
// semiring.
//
// This algorithm follows the specification given in the GAP Benchmark Suite:
// https://arxiv.org/abs/1508.03619

// For fastest results, the input matrix A should be stored in
// GxB_BY_COL format.  All entries in A must be equal to 1.

#include "LAGraph.h"

#define LAGRAPH_FREE_ALL {       \
    GrB_free(&d_out);            \
    GrB_free(&importance_vec);   \
    GrB_free(&pr);               \
    GrB_free(&oldpr);              \
    GrB_free (&op_diff) ;        \
};


void ddiff (void *z, const void *x, const void *y)
{
    float delta = (* ((float *) x)) - (* ((float *) y)) ;
    (*((float *) z)) = (delta > 0 ? delta : -delta);
}

GrB_Info LAGraph_pagerank3a // PageRank definition
(
 GrB_Vector *result,    // output: array of LAGraph_PageRank structs
 GrB_Matrix A,          // binary input graph, not modified
 float damping_factor, // damping factor
 unsigned long itermax, // maximum number of iterations
 int* iters             // output: number of iterations taken
 )
{
    GrB_Info info;
    GrB_Index n;
    GrB_Vector d_out = NULL ;
    GrB_Vector importance_vec = NULL ;
    GrB_Vector pr = NULL;
    GrB_BinaryOp op_diff = NULL ;
    GrB_Vector oldpr = NULL;

    LAGRAPH_OK(GrB_Matrix_nrows(&n, A));
    GrB_Index nvals;
    LAGRAPH_OK(GrB_Matrix_nvals(&nvals, A));

    // Matrix A row sum
    // Stores the outbound degrees of all vertices
    LAGRAPH_OK(GrB_Vector_new(&d_out, GrB_FP32, n));
    LAGRAPH_OK(GrB_reduce( d_out, NULL, NULL, GxB_PLUS_FP32_MONOID, A, NULL ));

    GrB_Index non_dangling ;
    LAGRAPH_OK (GrB_Vector_nvals (&non_dangling, d_out)) ;
    if (non_dangling < n)
    {
        LAGRAPH_ERROR ("Matrix has dangling nodes!", GrB_INVALID_VALUE) ;
    }

    // Iteration
    // Initialize PR vector
    LAGRAPH_OK(GrB_Vector_new(&pr, GrB_FP32, n));
    // Fill result vector with initial value (1 / |V|)
    LAGRAPH_OK(GrB_assign( pr, NULL, NULL, 1.0 / n, GrB_ALL, n, NULL ));

    LAGRAPH_OK(GrB_Vector_new(&importance_vec, GrB_FP32, n));

    // Teleport value
    const float teleport = (1 - damping_factor) / n;

    // create operator
    LAGRAPH_OK (GrB_BinaryOp_new (&op_diff, ddiff,  GrB_FP32, GrB_FP32, 
                GrB_FP32)) ;

    float  tol = 1e-4;
    float  rdiff = 1 ;       // so first iteration is always done

    for ((*iters) = 0 ; (*iters) < itermax && rdiff > tol ; (*iters)++)
    {
        // oldpr = pr; deep copy
        GrB_Vector_dup(&oldpr, pr);
 
        //
        // Importance calculation
        //

        // Divide previous PageRank with number of outbound edges
        // importance_vec = Pr/|d_out|
        LAGRAPH_OK(GrB_eWiseMult(  importance_vec, NULL, NULL, GrB_DIV_FP32,
                    pr, d_out, NULL ));

        // Multiply importance by damping factor, 
        // importance_vec = damping_factor . Pr/|d_out|
        // importance_vec .= damping_factor
        LAGRAPH_OK(GrB_assign( importance_vec, NULL, GrB_TIMES_FP32,
                    damping_factor, GrB_ALL, n, NULL ));

        //printf ("before mxv\n");
        //GxB_print(importance_vec, 3);
        
        // Calculate total PR of all inbound vertices
        // importance_vec = A' * importance_vec
        LAGRAPH_OK(GrB_mxv( importance_vec, NULL, NULL, GxB_PLUS_SECOND_FP32,
                    A, importance_vec, LAGraph_desc_tooo)) ;

        //GxB_print(importance_vec, 3);

        // PageRank summarization
        // Add teleport, importance_vec, and dangling_vec components together
        // pr = (1-df)/n
        LAGRAPH_OK(GrB_assign( pr, NULL, NULL, teleport, GrB_ALL, n, NULL ));
        // pr += importance_vec
        LAGRAPH_OK(GrB_eWiseAdd( pr, NULL, NULL, GxB_PLUS_FP32_MONOID,
                    pr, importance_vec, NULL ));

        //printf ("after mxv\n");
        //GxB_print(pr, 3);
        //----------------------------------------------------------------------
        // rdiff = sum ((pr-oldpr).^2)
        //----------------------------------------------------------------------
        // oldpr = abs(pr-oldpr)
        LAGRAPH_OK (GrB_eWiseAdd (oldpr, NULL, NULL, op_diff, oldpr, pr, NULL));
        // ridiff = sum (oldpr)
        LAGRAPH_OK (GrB_reduce (&rdiff, NULL, 
                    GxB_PLUS_FP32_MONOID, oldpr, NULL)) ;
        //printf("iters %d  rdiff=%f\n",*iters, rdiff);
        GrB_free (&oldpr) ;
   }

    (*result) = pr;
    pr = NULL ;
    LAGRAPH_FREE_ALL ;
    return (GrB_SUCCESS);
}
