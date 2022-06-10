//------------------------------------------------------------------------------
// LAGr_SquareClustering: vertex square-clustering
//------------------------------------------------------------------------------

// LAGraph, (c) 2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Erik Welch, NVIDIA.

//------------------------------------------------------------------------------

// TODO: summarize.

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
    GrB_free (&Q) ;                 \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
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

    // out_degrees assigned to diagonal matrix
    GrB_Matrix D = NULL ;

    // plus_pair(A @ A.T).new(mask=~D.S)
    GrB_Matrix P2 = NULL ;

    // Temporary workspace matrix (int64)
    GrB_Matrix Q = NULL ;

    GrB_Vector d_out = G->out_degree ;
    GrB_Matrix A = G->A ;
    GrB_Index n = 0 ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_ASSERT (square_clustering != NULL, GrB_NULL_POINTER) ;
    (*square_clustering) = NULL ;

    LG_ASSERT_MSG (d_out != NULL,
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
    GRB_TRY (GrB_Matrix_diag (&D, d_out, 0)) ;

    // We'll use `P2 = plus_pair(A @ A.T).new(mask=~D.S)` throughout.
    // We use ~D.S as a mask so P2 won't have values along the diagonal.
    GRB_TRY (GrB_Matrix_new (&P2, GrB_INT64, n, n)) ;
    GRB_TRY (GrB_mxm (P2, D, NULL, GxB_PLUS_PAIR_INT64, A, A, GrB_DESC_SCT1)) ;

    // Now compute the number of squares (the numerator).  We cound squares
    // based on https://arxiv.org/pdf/2007.11111.pdf (sigma_12, c_4).
    //     Q = P2 * (P2 - 1)
    //     squares = Q.reduce_rowwise() / 2  (and drop zeros)
    GRB_TRY (GrB_Matrix_new (&Q, GrB_INT64, n, n)) ;
    GRB_TRY (GrB_Matrix_apply_BinaryOp2nd_INT64 (Q, NULL, NULL, GrB_MINUS_INT64,
        P2, 1, NULL)) ;
    GRB_TRY (GrB_Matrix_eWiseMult_BinaryOp (Q, NULL, NULL, GrB_TIMES_INT64, Q,
        P2, NULL)) ;
    GRB_TRY (GrB_Vector_new (&squares, GrB_INT64, n)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (squares, NULL, NULL,
        GrB_PLUS_MONOID_INT64, Q, NULL)) ;
    // Divide by 2, and use squares as value mask to drop zeros
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64 (squares, squares, NULL,
        GrB_DIV_INT64, squares, 2, GrB_DESC_R));

    // Denominator is thought of as total number of squares that could exist.
    // We use the definition from https://arxiv.org/pdf/0710.0117v1.pdf.
    // First three contributions will become negative in the final step.
    //
    // (1) Subtract 1 for each u and 1 for each w for all combos:
    //     denom = d_out * (d_out - 1)
    GRB_TRY (GrB_Vector_new (&denom, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_INT64 (denom, NULL, NULL,
        GrB_MINUS_INT64, d_out, 1, NULL)) ;
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (denom, NULL, NULL, GrB_TIMES_INT64,
        denom, d_out, NULL)) ;

    // (2) Subtract the number of squares (will become negative, so add here):
    //     denom = denom + squares
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (denom, NULL, NULL, GrB_PLUS_INT64,
        denom, squares, NULL)) ;

    // (3) Subtract 1 for each edge where u-w or w-u are connected.
    // In other words, triangles.  Use P2, since we already have it.
    //     Q = first(P2 & A)
    //     denom += Q.reduce_rowwise()
    GRB_TRY (GrB_Matrix_eWiseMult_BinaryOp (Q, NULL, NULL, GrB_FIRST_INT64, P2,
        A, NULL)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (denom, NULL, GrB_PLUS_INT64,
        GrB_PLUS_MONOID_INT64, Q, NULL)) ;

    // The main contribution to the denominator:
    //     degrees[u] + degrees[w] for each u-w combo.
    // This is the only positive term.
    // We subtract all other terms from this one, hence rminus.
    //     Q = plus_pair(A @ P2.T).new(mask=A.S)
    //     Q = any_times(Q @ D)
    //     denom(rminus) = Q.reduce_rowwise()
    GRB_TRY (GrB_mxm (Q, A, NULL, GxB_PLUS_PAIR_INT64, A, P2, GrB_DESC_RST1)) ;
    GRB_TRY (GrB_mxm (Q, NULL, NULL, GxB_ANY_TIMES_INT64, Q, D, NULL)) ;
    GRB_TRY (GrB_Matrix_reduce_Monoid (denom, NULL, GxB_RMINUS_INT64,
        GrB_PLUS_MONOID_INT64, Q, NULL)) ;

    // Almost done!  Now compute the final result:
    //     square_clustering = squares / denom
    GRB_TRY (GrB_Vector_new (&r, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_eWiseMult_BinaryOp (r, NULL, NULL, GrB_DIV_FP64,
        squares, denom, NULL)) ;

    (*square_clustering) = r ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
