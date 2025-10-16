//------------------------------------------------------------------------------
// LAGraph_Jaccard - parallel jaccard similarity
// our second approach to compute jaccard similarity
//------------------------------------------------------------------------------

#define LG_FREE_WORK                           \
{                                              \
    GrB_free(&deg);                             \
    GrB_free(&Au);                              \
    GrB_free(&R);                               \
    GrB_free(&J);                               \
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
    // input: not modified
    LAGraph_Graph G,
	bool all_pairs, 
    char *msg
)
{
    GrB_Matrix B = NULL, J = NULL, Au = NULL, R = NULL;    
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

    GrB_Semiring semiring = NULL;

    GrB_Matrix A = G->A ;
    struct timeval stop, start, st_deg, en_deg, st_intersection, en_intersection, st_union, en_union, en_select;

    gettimeofday(&start, NULL);

	GRB_TRY( GrB_Matrix_nrows (&n, A) );
    //--------------------------------------------------------------------------
    // calculating out degree matrix deg
    //--------------------------------------------------------------------------
    gettimeofday(&st_deg, NULL);
    GrB_Type int_type  = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    GRB_TRY (GrB_Vector_new (&deg, int_type, n)) ;
	GRB_TRY( GrB_reduce(deg, NULL, NULL, (int_type == GrB_INT64) ? GrB_PLUS_INT64 : GrB_PLUS_INT32, A, NULL));
    gettimeofday(&en_deg, NULL);

	// GRB_TRY (LAGraph_Vector_Print (deg, LAGraph_COMPLETE_VERBOSE, stdout, msg)) ;
    //--------------------------------------------------------------------------
    // B is intersection matrix 
    //--------------------------------------------------------------------------
    gettimeofday(&st_intersection, NULL);
	GRB_TRY(GrB_Matrix_new(&B, GrB_FP32, n, n));

	//make a copy of A
	GRB_TRY(GrB_Matrix_new(&Au, GrB_UINT32, n, n));
	GRB_TRY( GrB_select(Au, NULL, NULL, GrB_VALUENE_BOOL, A, 0, NULL));
	GRB_TRY(GrB_mxm(B, all_pairs ? NULL : A, NULL, GxB_PLUS_TIMES_UINT32, Au, Au, NULL));
	
	GRB_TRY( GrB_select(B, NULL, NULL, GrB_TRIU, B, (int64_t)1, NULL));
    GRB_TRY (GrB_Matrix_wait (B, GrB_COMPLETE)) ;
    gettimeofday(&en_select, NULL);

	
    //--------------------------------------------------------------------------
    // S is jaccard index
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
    gettimeofday(&en_union, NULL);
	GRB_TRY (GrB_Matrix_wait (B, GrB_COMPLETE)) ;

    gettimeofday(&stop, NULL);

    long duration_usec = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
    printf("execute time is %lu microseconds\n", duration_usec);
    printf("degree calculation time is %lu microseconds\n", ((en_deg.tv_sec - st_deg.tv_sec) * 1000000 + en_deg.tv_usec - st_deg.tv_usec));
	printf("intesection and select time is %lu microseconds\n", 
    ((en_select.tv_sec - st_intersection.tv_sec) * 1000000 + en_select.tv_usec - st_intersection.tv_usec));
    printf("union time is %lu microseconds\n", ((en_union.tv_sec - st_union.tv_sec) * 1000000 + en_union.tv_usec - st_union.tv_usec));
    
   
	// GRB_TRY (LAGraph_Matrix_Print (B, LAGraph_COMPLETE, stdout, msg)) ;
    (*JC) = B;
	B = NULL;           
	LG_FREE_WORK;   

    return (GrB_SUCCESS) ;
}
