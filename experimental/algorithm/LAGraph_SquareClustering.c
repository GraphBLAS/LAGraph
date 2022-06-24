//------------------------------------------------------------------------------
// LAGraph_SquareClustering: vertex square-clustering
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Erik Welch, NVIDIA.

//------------------------------------------------------------------------------

// TODO: summarize.  This calculates `P2 = A @ A.T`, which may be very large.

// https://networkx.org/documentation/stable/reference/algorithms/generated/\
//         networkx.algorithms.cluster.square_clustering.html
// https://arxiv.org/pdf/2007.11111.pdf
// https://arxiv.org/pdf/0710.0117v1.pdf

//------------------------------------------------------------------------------

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&squares) ;           \
    GrB_free (&denom) ;             \
    GrB_free (&D) ;                 \
    GrB_free (&P2) ;                \
    GrB_free (&uw_degrees) ;        \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (&Tri) ;               \
    GrB_free (&r) ;                 \
}

#include <LAGraph.h>
#include <LAGraphX.h>
#include <LG_internal.h>  // from src/utility

int LAGraph_SquareClustering
(
    // outputs:
    GrB_Vector *square_clustering,
    // inputs:
    LAGraph_Graph G,
    char *msg
)
{
    LG_CLEAR_MSG ;

    // The number of squares each node is part of
    GrB_Vector squares = NULL ;

    // Thought of as the total number of possible squares for each node
    GrB_Vector denom = NULL ;

    // Final result: the square coefficients for each node (squares / denom)
    GrB_Vector r = NULL ;

    // Used by pure GrB version; the largest factor of the denominator
    GrB_Vector uw_degrees = NULL ;

    // out_degrees assigned to diagonal matrix
    GrB_Matrix D = NULL ;

    // P2 = plus_pair(A @ A.T).new(mask=~D.S)
    // Then used as a temporary workspace matrix (int64)
    GrB_Matrix P2 = NULL ;

    // Triangles: first(P2 & A)
    GrB_Matrix Tri = NULL ;

    GrB_Vector deg = G->out_degree ;
    GrB_Matrix A = G->A ;
    GrB_Index n = 0 ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_ASSERT (square_clustering != NULL, GrB_NULL_POINTER) ;
    (*square_clustering) = NULL ;

    LG_ASSERT_MSG (deg != NULL,
        LAGRAPH_NOT_CACHED, "G->out_degree is required") ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    LG_ASSERT_MSG ((G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->is_symmetric_structure == LAGraph_TRUE)),
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be known to be symmetric") ;

    // # of nodes
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    // out_degrees as a diagonal matrix.  We use this twice:
    // 1) as a mask to ignore the diagonal elements when computing `A @ A.T`
    // 2) and to multiply each column of a matrix by the degrees.
    #if LAGRAPH_SUITESPARSE
        #if GxB_IMPLEMENTATION >= GxB_VERSION (7,0,0)
        // SuiteSparse 7.x and later:
        GRB_TRY (GrB_Matrix_diag(&D, deg, 0)) ;
        #else
        // SuiteSparse 6.x and earlier, which had the incorrect signature:
        GRB_TRY (GrB_Matrix_new(&D, GrB_INT64, n, n)) ;
        GRB_TRY (GrB_Matrix_diag(D, deg, 0)) ;
        #endif
    #else
    // standard GrB:
    GRB_TRY (GrB_Matrix_diag(&D, deg, 0)) ;
    #endif

    // We'll use `P2 = plus_pair(A @ A.T).new(mask=~D.S)` throughout.
    // We use ~D.S as a mask so P2 won't have values along the diagonal.
    GRB_TRY (GrB_Matrix_new (&P2, GrB_INT64, n, n)) ;
    GRB_TRY (GrB_mxm (P2, D, NULL, LAGraph_plus_one_int64, A, A, GrB_DESC_SCT1)) ;

    // Denominator is thought of as total number of squares that could exist.
    // We use the definition from https://arxiv.org/pdf/0710.0117v1.pdf.
    // First three contributions will become negative in the final step.
    //
    // (1) Subtract 1 for each u and 1 for each w for all combos:
    //     denom = deg * (deg - 1)
    GRB_TRY (GrB_Vector_new (&denom, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64 (denom, NULL, NULL,
        GrB_MINUS_INT64, deg, 1, NULL)) ;
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (denom, NULL, NULL, GrB_TIMES_INT64,
        denom, deg, NULL)) ;

    // (2) Subtract 1 for each edge where u-w or w-u are connected.
    // In other words, triangles.  Use P2, since we already have it.
    //     Tri = first(P2 & A)
    //     denom += Tri.reduce_rowwise()
    GRB_TRY (GrB_Matrix_new (&Tri, GrB_INT64, n, n)) ;
    GRB_TRY (GrB_Matrix_eWiseMult_BinaryOp (Tri, NULL, NULL, GrB_FIRST_INT64, P2,
        A, NULL)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (denom, NULL, GrB_PLUS_INT64,
        GrB_PLUS_MONOID_INT64, Tri, NULL)) ;
    GrB_free (&Tri) ;

    // Now compute the number of squares (the numerator).  We count squares
    // based on https://arxiv.org/pdf/2007.11111.pdf (sigma_12, c_4).
    //     P2 *= P2 - 1
    //     squares = P2.reduce_rowwise() / 2  (and drop zeros)
    GRB_TRY (GrB_Matrix_apply_BinaryOp2nd_INT64 (P2, NULL, GrB_TIMES_INT64,
        GrB_MINUS_INT64, P2, 1, NULL)) ;
    GRB_TRY (GrB_Vector_new (&squares, GrB_INT64, n)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (squares, NULL, NULL,
        GrB_PLUS_MONOID_INT64, P2, NULL)) ;
    // Divide by 2, and use squares as value mask to drop zeros
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64 (squares, squares, NULL,
        GrB_DIV_INT64, squares, 2, GrB_DESC_R)) ;

    // (3) Subtract the number of squares (will become negative, so add here):
    //     denom = denom + squares
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (denom, NULL, NULL, GrB_PLUS_INT64,
        denom, squares, NULL)) ;

    // The main contribution to the denominator:
    //     degrees[u] + degrees[w] for each u-w combo.
    // This is the only positive term.
    // We subtract all other terms from this one, hence rminus.
    //     P2 = plus_pair(A @ P2.T).new(mask=A.S)
    //     P2 = any_times(P2 @ D)
    //     denom(rminus) = P2.reduce_rowwise()
    GRB_TRY (GrB_mxm (P2, A, NULL, LAGraph_plus_one_int64, A, P2,
        GrB_DESC_RST1)) ;
    GRB_TRY (GrB_mxm (P2, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, P2, D,
        NULL)) ;
    #if LAGRAPH_SUITESPARSE
    GRB_TRY (GrB_Matrix_reduce_Monoid (denom, NULL, GxB_RMINUS_INT64,
        GrB_PLUS_MONOID_INT64, P2, NULL)) ;
    #else
    GRB_TRY (GrB_Vector_new (&uw_degrees, GrB_INT64, n)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (uw_degrees, NULL, NULL,
        GrB_PLUS_MONOID_INT64, P2, NULL)) ;
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (denom, NULL, NULL,
        GrB_MINUS_INT64, uw_degrees, denom, NULL)) ;
    #endif

    // Almost done!  Now compute the final result:
    //     square_clustering = squares / denom
    GRB_TRY (GrB_Vector_new (&r, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (r, NULL, NULL, GrB_DIV_FP64,
        squares, denom, NULL)) ;

    (*square_clustering) = r ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
