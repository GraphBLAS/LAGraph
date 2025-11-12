//------------------------------------------------------------------------------
// LAGraph_Jaccard - parallel jaccard similarity
//------------------------------------------------------------------------------

#define LG_FREE_WORK                           \
{                                              \
    GrB_free(&deg);                             \
    GrB_free(&Au);                              \
    GrB_free(&R);                               \
    GrB_free(&B);                               \
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
#include <LAGraph.h>
#include "LAGraphX.h" 
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int LAGraph_Jaccard // a simple algorithm, just for illustration
(
    // output
    GrB_Matrix *JC,
    // input:
    LAGraph_Graph G,
	bool all_pairs, 
    char *msg
)
{
    GrB_Matrix B = NULL, Au = NULL, R = NULL;    
    GrB_Index n;
    GrB_Vector deg = NULL;
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    LG_CLEAR_MSG ;

	int nthreads = omp_get_max_threads() ;
	printf("num of threads %d\n", nthreads);

    LG_ASSERT (JC != NULL, GrB_NULL_POINTER) ;
    (*JC) = NULL ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    
    GrB_Matrix A = G->A ;
	GRB_TRY( GrB_Matrix_nrows (&n, A) );
    //--------------------------------------------------------------------------
    // calculating degree vector deg
    //--------------------------------------------------------------------------
    
    GrB_Type int_type  = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    GRB_TRY (GrB_Vector_new (&deg, int_type, n)) ;
	GRB_TRY( GrB_reduce(deg, NULL, NULL, (int_type == GrB_INT64) ? GrB_PLUS_INT64 : GrB_PLUS_INT32, A, NULL));
    
	// GRB_TRY (LAGraph_Vector_Print (deg, LAGraph_COMPLETE_VERBOSE, stdout, msg)) ;
    //--------------------------------------------------------------------------
    // B is intersection matrix 
    //--------------------------------------------------------------------------
    
	GRB_TRY(GrB_Matrix_new(&B, GrB_FP32, n, n));

	//make a copy of A
	GRB_TRY(GrB_Matrix_new(&Au, GrB_UINT32, n, n));
	GRB_TRY( GrB_select(Au, NULL, NULL, GrB_VALUENE_BOOL, A, 0, NULL));
	GRB_TRY(GrB_mxm(B, all_pairs ? NULL : A, NULL, GxB_PLUS_TIMES_UINT32, Au, Au, NULL));
	
	GRB_TRY( GrB_select(B, NULL, NULL, GrB_TRIU, B, (int64_t)1, NULL));
    GRB_TRY (GrB_Matrix_wait (B, GrB_COMPLETE)) ;
    gettimeofday(&en_select, NULL);

	
    //--------------------------------------------------------------------------
    // R has summation of degree of corresponding nodes
    // B is jaccard index B <- B / (R-B)
    //--------------------------------------------------------------------------
	// assign deg //
	GrB_Matrix_new(&R, GrB_FP32, n, n);
	for (GrB_Index j = 0; j < n; j++) {
		GRB_TRY( GrB_assign(R, NULL, NULL, deg, (GrB_Index*) GrB_ALL, n, j, NULL));
	}
	// assign deg into every row and ADD to current R: R(i, :) += v^T
	for (GrB_Index i = 0; i < n; i++) {
		GRB_TRY( GrB_assign(R, NULL, GrB_PLUS_INT32, deg, i, (GrB_Index*) GrB_ALL, n, GrB_DESC_T0));
	}

	GRB_TRY(  GrB_eWiseAdd(R, B, NULL, GrB_MINUS_FP32, R, B, NULL) );
	GRB_TRY(  GrB_eWiseMult(B, NULL, NULL, GrB_DIV_FP32, B, R, NULL) );
	GRB_TRY (GrB_Matrix_wait (B, GrB_COMPLETE)) ;
   
	// GRB_TRY (LAGraph_Matrix_Print (B, LAGraph_COMPLETE, stdout, msg)) ;
    (*JC) = B;
	B = NULL;           
	LG_FREE_WORK;   

    return (GrB_SUCCESS) ;
}
