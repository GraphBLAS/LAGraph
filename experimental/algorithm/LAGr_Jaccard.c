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

// TODO: add a description, and citations

// References:
// (1) (your paper)
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

    // TODO: what if deg is sparse?  We could do
    #if 0
    if nvals(deg) < n
        t = deg (via GrB_Vector_dup)
        t < !t, struct> = 0
        d = t
    else
        d = deg alias
    when done, free the t vector
    #endif

    // D is degree matrix
    GRB_TRY(GrB_Matrix_diag(&D, deg, 0));   // use d here
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
