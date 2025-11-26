//------------------------------------------------------------------------------
// LAGraph_Jaccard - parallel jaccard similarity
//------------------------------------------------------------------------------

#define LG_FREE_WORK                           \
{                                              \
    GrB_free(&R);                               \
    GrB_free(&B);                               \
	GrB_free(&D);                               \
}

#define LG_FREE_ALL                            \
{                                              \
    LG_FREE_WORK ;                             \
    if (JC != NULL && *JC != NULL)             \
    {                                          \
        GrB_free(JC);                           \
    }                                          \
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
    GrB_Matrix B = NULL, R = NULL, D = NULL;    
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
	GRB_TRY(GrB_Matrix_new(&B, GrB_FP64, n, n));

	if (!all_pairs)
	{
		GRB_TRY( GrB_select(B, NULL, NULL, GrB_TRIU, A, (int64_t)0, NULL));
	}

	GRB_TRY(GrB_mxm(B, all_pairs ? NULL : B, NULL, LAGraph_plus_one_uint32, A, A, GrB_DESC_S));
	
	if (all_pairs)
	{
		GRB_TRY( GrB_select(B, NULL, NULL, GrB_TRIU, B, (int64_t)0, NULL));
	}
	
    //--------------------------------------------------------------------------
    // R has summation of degree of corresponding nodes
    // B is jaccard index B <- B / (R-B)
    //--------------------------------------------------------------------------
	// D is degree matrix
	GRB_TRY(GrB_Matrix_diag(&D, deg, 0));
	GrB_Matrix_new(&R, int_type, n, n);

	// R = B*D  -> R_ij = deg_j - b_ij
	GRB_TRY(GrB_mxm(R, all_pairs ? NULL : B, NULL, (int_type == GrB_INT64) ? GxB_PLUS_RMINUS_INT64 : GxB_PLUS_RMINUS_INT32, B, D, GrB_DESC_S));
	// R = B*D  -> R_ij = deg_j - b_ij	
	GRB_TRY(GrB_mxm(R, all_pairs ? NULL : B, NULL, (int_type == GrB_INT64) ? GxB_PLUS_PLUS_INT64 : GxB_PLUS_PLUS_INT32, D, R, GrB_DESC_S));	

	GRB_TRY(  GrB_eWiseMult(B, NULL, NULL, (int_type == GrB_INT64) ? GrB_DIV_FP64 : GrB_DIV_FP32, B, R, NULL) );
    (*JC) = B;
	B = NULL;           
	LG_FREE_WORK;   
	
    return (GrB_SUCCESS) ;
}
