//------------------------------------------------------------------------------
// LAGraph_Jaccard - parallel jaccard similarity
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2025 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Elaheh Hassani and Tim Davis, Texas A&M University

//------------------------------------------------------------------------------
// LAGr_Jaccard: compute Jaccard similarity (weight) coefficients for an undirected graph.
//
// Inputs:
// G             - a valid LAGraph_Graph with:
//                    * G->A structurally symmetric (undirected)
//                    * no self-edges
//                    * G->out_degree cached (no explicit zeros)
// all_pairs     - if true, compute Jaccard for all pairs of nodes in the graph;
//                 if false, compute only for neighboring nodes (current edges).
//
// Output:
// JC            - an n-by-n matrix of Jaccard coefficients. JC is always
//                 returned as an upper-triangular matrix. When all_pairs is
//                 true, the diagonal entries are equal to 1. When all_pairs is
//                 false, the diagonal is structurally empty.
//
// Jaccard similarity measures the overlap of neighbor sets.
// For a pair of nodes i and j in the graph G:
//         JC(i,j) = |N(i) ∩ N(j)| / |N(i) ∪ N(j)|
// where N(i) is the set of neighbors of i in G.
//
// The computation proceeds in two stages.  First, an “intersection matrix” B
// is formed where B(i,j) = |N(i) ∩ N(j)|.  
// Second, a degree-based denominator matrix R is formed where
// R(i,j) = deg(i) + deg(j) − B(i,j).
// Lastly, the output matrix JC is computed as the element-wise ratio
// JC = B ./ R.
//

//------------------------------------------------------------------------------
// References:
// (1) "Parallel Algorithms for Computing Jaccard Weights on Graphs using Linear Algebra," 
//     in Proc. IEEE High Performance Extreme Computing Conference (HPEC), 2023.
//     https://doi.org/10.1109/HPEC58863.2023.10363558
//
// (2) https://en.wikipedia.org/wiki/Jaccard_index


#define LG_FREE_WORK                           \
{                                              \
    GrB_free(&R);                              \
    GrB_free(&B);                              \
    GrB_free(&D);                              \
    GrB_free (&M) ;                            \
}

#define LG_FREE_ALL                            \
{                                              \
    LG_FREE_WORK ;                             \
}

#include "LG_internal.h"

int LAGr_Jaccard
(
    // output
    GrB_Matrix *JC,
    // input:
    LAGraph_Graph G,
    bool all_pairs, 
    char *msg
)
{
    GrB_Matrix B = NULL, R = NULL, D = NULL, M = NULL ;    
    GrB_Index n;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    LG_CLEAR_MSG ;
    // error if G is directed (OK if G is directed but G->A is symmetric in structure
    // error if G has self edges, or unknown
    LG_ASSERT (JC != NULL, GrB_NULL_POINTER) ;
    (*JC) = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, "G->nself_edges must be zero") ;
    LG_ASSERT_MSG ((G->kind == LAGraph_ADJACENCY_UNDIRECTED || (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->is_symmetric_structure == LAGraph_TRUE)),
        LAGRAPH_SYMMETRIC_STRUCTURE_REQUIRED,
        "G->A must be known to be symmetric") ;
    //--------------------------------------------------------------------------
    // degree vector deg
    //--------------------------------------------------------------------------
    GrB_Vector deg = G->out_degree ;
    LG_ASSERT_MSG (deg != NULL, LAGRAPH_NOT_CACHED, "G->out_degree is required") ;
    GrB_Matrix A = G->A ;
    GRB_TRY( GrB_Matrix_nrows (&n, A) );
    GrB_Type int_type  = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    
    //--------------------------------------------------------------------------
    // B is intersection matrix 
    //--------------------------------------------------------------------------

    // B(i,j) is the size of the intersection of the pattern of A(i,:) and
    // A(:,j).  If all_pairs is true, B is computed for all entries in A^2.
    // Otherwise, it is computed just for entries in triu(A).

    // The final output matrix (B) is always upper triangular, in both cases.

    GRB_TRY(GrB_Matrix_new(&B, GrB_FP64, n, n));

    if (all_pairs)
    {
        // B = triu (A*A)
        GRB_TRY(GrB_mxm(B, NULL, NULL, LAGraph_plus_one_fp64, A, A, NULL));
        GRB_TRY( GrB_select(B, NULL, NULL, GrB_TRIU, B, (int64_t)0, NULL));
    }
    else
    {
        // B<triu(A)> = A*A'
        GRB_TRY(GrB_Matrix_new(&M, GrB_BOOL, n, n));
        GRB_TRY( GrB_select(M, NULL, NULL, GrB_TRIU, A, (int64_t)0, NULL));
        GRB_TRY(GrB_mxm(B, M, NULL, LAGraph_plus_one_fp64, A, A, GrB_DESC_ST1));
        GrB_free (&M) ;
    }

    //--------------------------------------------------------------------------
    // R has summation of degree of corresponding nodes
    // B is jaccard index B <- B / (R-B)
    //--------------------------------------------------------------------------

    // If deg vectors is sparse, make it dense
	GrB_Index deg_nnz = 0;
	GrB_Vector d = NULL;
    GRB_TRY(GrB_Vector_nvals(&deg_nnz, deg));
	if (deg_nnz < n ) {
		GrB_Vector t = NULL;
		GRB_TRY (GrB_Vector_dup(&t, deg));
		GRB_TRY (GrB_assign (t, t, NULL, 0, GrB_ALL, n, GrB_DESC_SC)) ;
		GRB_TRY (GrB_Vector_dup(&d, t));
	}
	else {
		GRB_TRY (GrB_Vector_dup(&d, deg));
	}

    // D is degree matrix
    GRB_TRY(GrB_Matrix_diag(&D, d, 0));   // use d here
    GrB_Matrix_new(&R, int_type, n, n);

    // R = B*D  -> R_ij = deg_j - b_ij
    GRB_TRY(GrB_mxm(R, NULL, NULL, (int_type == GrB_INT64) ? GxB_PLUS_RMINUS_INT64 : GxB_PLUS_RMINUS_INT32, B, D, GrB_DESC_S));
    // R = D*R  -> R_ij = deg_i + r_ij    
    GRB_TRY(GrB_mxm(R, NULL, NULL, (int_type == GrB_INT64) ? GxB_PLUS_PLUS_INT64 : GxB_PLUS_PLUS_INT32, D, R, GrB_DESC_S));    

    GRB_TRY(  GrB_eWiseMult(B, NULL, NULL, GrB_DIV_FP64, B, R, NULL) );
    (*JC) = B;
    B = NULL;           
    LG_FREE_WORK;   
    
    return (GrB_SUCCESS) ;
}
