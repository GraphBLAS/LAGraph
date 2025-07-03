//------------------------------------------------------------------------------
// LAGr_MaxFlow: max flow
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Darin Peries and Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// FIXME: the GrB_ and GRB_ prefixes should be reserved for GraphBLAS itself,
// and not used for new matrices and operators here.  Rename them.

#include <LAGraphX.h>
#include "LG_internal.h"
#include <LAGraph.h>

#if LG_SUITESPARSE_GRAPHBLAS_V10

//------------------------------------------------------------------------------
// LG_augment_maxflow
//------------------------------------------------------------------------------

// LG_augment_maxflow is a function used to sum the current excess flow of the
// sink into the output variable f for each iteration.

#undef  LG_FREE_ALL
#define LG_FREE_ALL ;

static GrB_Info LG_augment_maxflow
(
    double *f,                  // total maxflow from src to sink
    GrB_Vector e,               // excess vector, of type double
    GrB_Index sink,             // sink node
    GrB_Vector src_and_sink,    // mask vector, with just [src sink]
    GrB_Index *n_active,        // # of active nodes
    char *msg
)
{
    // e_sink = e (sink)
    double e_sink = 0;
    GrB_Info info = GrB_Vector_extractElement(&e_sink, e, sink); //if value at sink
    GRB_TRY (info) ;
    if (info == GrB_SUCCESS)
    {
        // e(sink) is present
        (*f) += e_sink; 
    }

    // TODO: what if e is tiny?  Do we need a tol parameter,
    // and replace all "e > 0" and "r > 0" comparisons with
    // (e > tol), throughout the code?

    // e<![src,sink]> = select e where (e > 0)
    GRB_TRY(GrB_select(e, src_and_sink, NULL, GrB_VALUEGT_FP64, e, 0, GrB_DESC_RSC));   
    GRB_TRY(GrB_Vector_nvals(n_active, e));
}

//------------------------------------------------------------------------------

#undef LG_FREE_WORK
#undef LG_FREE_ALL

#define LG_FREE_WORK                                                           \
  {                                                                            \
    GrB_free(&GrB_FlowEdge);                                                   \
    GrB_free(&GrB_CompareTuple);                                               \
    GrB_free(&GrB_ResultTuple);                                                \
    GrB_free(&e);                                                              \
    GrB_free(&d);                                                              \
    GrB_free(&theta);                                                          \
    GrB_free(&R);                                                              \
    GrB_free(&delta);                                                          \
    GrB_free(&delta_vec);                                                      \
    GrB_free(&delta_mat);                                                      \
    GrB_free(&map);                                                            \
    GrB_free(&y);                                                              \
    GrB_free(&yd);                                                             \
    GrB_free(&src_and_sink);                                                   \
    GrB_free(&Jvec);                                                           \
    GrB_free(&GrB_Prune);                                                      \
    GrB_free(&GrB_UpdateFlows);                                                \
    GrB_free(&GrB_UpdateHeight);                                               \
    GrB_free(&GrB_extractFlows);                                               \
    GrB_free(&GrB_MxeIndexMult);                                               \
    GrB_free(&GrB_MxeMult);                                                    \
    GrB_free(&GrB_MxeAdd);                                                     \
    GrB_free(&GrB_MxeAddMonoid);                                               \
    GrB_free(&GrB_MxeSemiring);                                                \
    GrB_free(&GrB_extractJ);                                                   \
    GrB_free(&GrB_CreateCompareVec);                                           \
    GrB_free(&GrB_RxdSemiring);                                                \
    GrB_free(&GrB_RxdAdd);                                                     \
    GrB_free(&GrB_RxdAddMonoid);                                               \
    GrB_free(&GrB_RxdIndexMult);                                               \
    GrB_free(&GrB_RxdMult);                                                    \
    GrB_free(&GrB_InitForwardFlows);                                           \
    GrB_free(&GrB_InitBackwardFlows);                                          \
    GrB_free(&GrB_CreateResidualForward);                                      \
    GrB_free(&GrB_CreateResidualBackward);                                     \
    GrB_free(&zero);                                                           \
    GrB_free(&empty);                                                          \
    GrB_free(&Re);                                                             \
    GrB_free(&invariant);                                                      \
    GrB_free(&GrB_InvariantCheck);                                             \
    GrB_free(&check);                                                          \
    GrB_free(&GrB_extractYJ);                                                  \
    GrB_free(&extract_desc);                                                   \
    GrB_free(&residual_vec);						       \
    GrB_free(&GrB_MakeFlow);						       \
    GrB_free(&GrB_GetResidual);				                       \
  }


#define LG_FREE_ALL \
{ \
  LG_FREE_WORK; \
}

#define JIT_STR(f, var) char* var = #f; f

//casting for unary ops
#define F_UNARY(f) ((void (*)(void *, const void *))f)

// casting for index binary ops
#define F_INDEX_BINARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void *, GrB_Index, GrB_Index, const void *)) f)

// casting for binary op
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

// custom types
JIT_STR(typedef struct{
  double flow;
  double capacity;
  } MF_flowEdge;, GRB_FLOWEDGE_STR) 

JIT_STR(typedef struct{
  double residual;
  int64_t j;
  int64_t d;
  } MF_resultTuple64;, GRB_RESULTTUPLE_STR64)   

JIT_STR(typedef struct{
  double residual;
  int32_t j;
  int32_t d;
  } MF_resultTuple32;, GRB_RESULTTUPLE_STR32)   


JIT_STR(typedef struct{
  double residual;
  int64_t di;
  int64_t y_dmin;
  int64_t j;
  } MF_compareTuple64;, GRB_COMPARETUPLE_STR64) 

JIT_STR(typedef struct{
  double residual;
  int32_t di;
  int32_t y_dmin;
  int32_t j;
  int32_t unused;   /* to pad the struct to 24 bytes */
  } MF_compareTuple32;, GRB_COMPARETUPLE_STR32) // 24 bytes: padded


JIT_STR(void MF_CreateResidualForward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = (*y);
  }, GRB_CRF_STR)

JIT_STR(void MF_CreateResidualBackward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = 0;
  }, GRB_CRB_STR)


JIT_STR(void MF_RxdMult64(MF_resultTuple64 *z, const MF_flowEdge *x, GrB_Index ix,
			GrB_Index jx, const int64_t *y,
			GrB_Index iy, GrB_Index jy, const int64_t* theta) {
  double r = x->capacity - x->flow; 
  if(r > 0){
    z->d = *y;
    z->residual = r;
    z->j = jx;
  }
  else{
    z->d = INT64_MAX;
    z->residual = 0;
    z->j = -1;
  }
  }, GRB_RXDMULT_STR64)


JIT_STR(void MF_RxdMult32(MF_resultTuple32 *z, const MF_flowEdge *x, GrB_Index ix,
			GrB_Index jx, const int32_t *y,
			GrB_Index iy, GrB_Index jy, const int32_t* theta) {
  double r = x->capacity - x->flow; 
  if(r > 0){
    z->d = *y;
    z->residual = r;
    z->j = jx;
  }
  else{
    z->d = INT32_MAX;
    z->residual = 0;
    z->j = -1;
  }
}, GRB_RXDMULT_STR32)

JIT_STR(void MF_RxdAdd64(MF_resultTuple64 * z,
		       const MF_resultTuple64 * x, const MF_resultTuple64 * y) {
  if(x->d < y->d){
    (*z) = (*x) ;
  }
  else if(x->d > y->d){
    (*z) = (*y) ;
  }
  else{
    if(x->residual > y->residual){
      (*z) = (*x) ;
    }
    else if(x->residual < y->residual){
      (*z) = (*y) ;
    }
    else{
      if(x->j > y->j){
	(*z) = (*x);
      }
      else{
	(*z) = (*y) ;
      }
    }
  }
  }, GRB_RXDADD_STR64)

JIT_STR(void MF_RxdAdd32(MF_resultTuple32 * z,
		       const MF_resultTuple32 * x, const MF_resultTuple32 * y) {
  if(x->d < y->d){
    (*z) = (*x) ;
  }
  else if(x->d > y->d){
    (*z) = (*y) ;
  }
  else{
    if(x->residual > y->residual){
      (*z) = (*x) ;
    }
    else if(x->residual < y->residual){
      (*z) = (*y) ;
    }
    else{
      if(x->j > y->j){
	(*z) = (*x);
      }
      else{
	(*z) = (*y) ;
      }
    }
  }
  }, GRB_RXDADD_STR32)



JIT_STR(void MF_extractFlow64(double *z, const MF_resultTuple64 *x)
	{ (*z) = x->residual; }, GRB_EXTRACTFLOW_STR64)

JIT_STR(void MF_extractFlow32(double *z, const MF_resultTuple32 *x)
	{ (*z) = x->residual; }, GRB_EXTRACTFLOW_STR32)


JIT_STR(void MF_updateFlow(MF_flowEdge *z,
			   const MF_flowEdge *x, const double *y) {
  z->capacity = x->capacity;
  z->flow = x->flow + (*y); 
  }, GRB_UPDATEFLOWS_STR)


JIT_STR(void MF_updateHeight64(int64_t *z,
			     const int64_t *x, const MF_resultTuple64 *y) {
  if((*x) < y->d+1){
    (*z) = y->d + 1;
  }
  else {
    (*z) = (*x);
  }
  }, GRB_UPDATEHEIGHT_STR64)

JIT_STR(void MF_updateHeight32(int32_t *z,
			     const int32_t *x, const MF_resultTuple32 *y) {
  if((*x) < y->d+1){
    (*z) = y->d + 1;
  }
  else {
    (*z) = (*x);
  }
  }, GRB_UPDATEHEIGHT_STR32)


JIT_STR(void MF_extractJ64(int64_t *z, const MF_compareTuple64 *x) { (*z) = x->j; }, GRB_EXTRACTJ_STR64)

JIT_STR(void MF_extractJ32(int32_t *z, const MF_compareTuple32 *x) { (*z) = x->j; }, GRB_EXTRACTJ_STR32) 


JIT_STR(void MF_extractYJ64(int64_t *z, const MF_resultTuple64 *x) {
  (*z) = x->j;
  }, GRB_EXTRACTYJ_STR64)

JIT_STR(void MF_extractYJ32(int32_t *z, const MF_resultTuple32 *x) {
  (*z) = x->j;
  }, GRB_EXTRACTYJ_STR32)


JIT_STR(void MF_initForwardFlows(MF_flowEdge * z,
				 const MF_flowEdge * x, const MF_flowEdge * y){
  z->flow = y->flow + x->flow;  
  z->capacity = x->capacity;
  }, GRB_INITFLOWF_STR)


JIT_STR(void MF_initBackwardFlows(MF_flowEdge * z,
				  const MF_flowEdge * x, const MF_flowEdge * y){
  z->flow = x->flow - y->flow;  
  z->capacity = x->capacity;
  }, GRB_INITFLOWB_STR)

JIT_STR(void MF_MxeMult64(MF_resultTuple64 * z, const MF_compareTuple64 * x,
			GrB_Index ix, GrB_Index jx,
			const double * y, GrB_Index iy,
			GrB_Index jy, const int64_t* theta){
  if(x->di == x->y_dmin && (*y) > 0){ 
    if(ix < jx){
      z->d = x->y_dmin;
      z->residual = x->residual;
      z->j = x->j;
    } 
    else{
      z->d = INT64_MAX;
      z->residual = 0;
      z->j = -1;
    }
  }
  else if(x->di == x->y_dmin - 1 && (*y) > 0){
    z->d = INT64_MAX;
    z->residual = 0;
    z->j = -1;
  }
  else if(x->di <= x->y_dmin-1 || x->di == x->y_dmin+1 || x->di == x->y_dmin){
    z->d = x->y_dmin;
    z->residual = x->residual;
    z->j = x->j;
  }
  else{ 
    z->d = INT64_MAX;
    z->residual = 0;
    z->j = -1;
  }
  }, GRB_MXEMULT_STR64)

JIT_STR(void MF_MxeMult32(MF_resultTuple32 * z, const MF_compareTuple32 * x,
			GrB_Index ix, GrB_Index jx,
			const double * y, GrB_Index iy,
			GrB_Index jy, const int32_t* theta){
  if(x->di == x->y_dmin && (*y) > 0){   
    if(ix < jx){
      z->d = x->y_dmin;
      z->residual = x->residual;
      z->j = x->j;
    } 
    else{
      z->d = INT32_MAX;
      z->residual = 0;
      z->j = -1;
    }
  }
  else if(x->di == x->y_dmin - 1 && (*y) > 0){ 
    z->d = INT32_MAX;
    z->residual = 0;
    z->j = -1;
  }
  else if(x->di <= x->y_dmin-1 || x->di == x->y_dmin+1 || x->di == x->y_dmin){
    z->d = x->y_dmin;
    z->residual = x->residual;
    z->j = x->j;
  }
  else{ 
    z->d = INT32_MAX;
    z->residual = 0;
    z->j = -1;
  }
  }, GRB_MXEMULT_STR32)


// FIXME: x will always be non-NULL in the test below:
JIT_STR(void MF_MxeAdd64(MF_resultTuple64 * z,
		       const MF_resultTuple64 * x, const MF_resultTuple64 * y){
  if(x != NULL){
    (*z) = (*y) ;
  }
  else{
    (*z) = (*x) ;
  }
  }, GRB_MXEADD_STR64)

// FIXME: x will always be non-NULL in the test below:
JIT_STR(void MF_MxeAdd32(MF_resultTuple32 * z,
		       const MF_resultTuple32 * x, const MF_resultTuple32 * y){
  if(x != NULL){
    (*z) = (*y) ;
  }
  else{
    (*z) = (*x) ;
  }
  }, GRB_MXEADD_STR32)


JIT_STR(void MF_CreateCompareVec64(MF_compareTuple64 *comp,
				 const MF_resultTuple64 *res, const int64_t *height) {
  comp->di = (*height);
  comp->j = res->j;
  comp->residual = res->residual;
  comp->y_dmin = res->d;
  }, GRB_CREATECOMPVEC_STR64)

JIT_STR(void MF_CreateCompareVec32(MF_compareTuple32 *comp,
				 const MF_resultTuple32 *res, const int32_t *height) {
  comp->di = (*height);
  comp->j = res->j;
  comp->residual = res->residual;
  comp->y_dmin = res->d;
  comp->unused = 0 ;
  }, GRB_CREATECOMPVEC_STR32)


JIT_STR(void MF_Prune64(bool * z, const MF_resultTuple64 * x,
		      GrB_Index ix, GrB_Index jx, const int64_t * theta){
  *z = (x->j != *theta) ;
  }, GRB_PRUNE_STR64)

JIT_STR(void MF_Prune32(bool * z, const MF_resultTuple32 * x,
		      GrB_Index ix, GrB_Index jx, const int32_t * theta){
  *z = (x->j != *theta) ;
  }, GRB_PRUNE_STR32)


JIT_STR(void MF_MakeFlow(MF_flowEdge * flow_edge, const double * flow){
  flow_edge->capacity = 0;
  flow_edge->flow = (*flow);
  }, GRB_MAKEF_STR)

#ifdef DBG
JIT_STR(void MF_CheckInvariant64(bool *z, const int64_t *height,
			       const MF_resultTuple64 *result) {
  (*z) = ((*height) == result->d+1);
  }, GRB_INV_STR64)

JIT_STR(void MF_CheckInvariant32(bool *z, const int32_t *height,
			       const MF_resultTuple32 *result) {
  (*z) = ((*height) == result->d+1);
  }, GRB_INV_STR32)
#endif
  
JIT_STR(void MF_getResidual(double * res, const MF_flowEdge * flow_edge){
    (*res) = flow_edge->capacity - flow_edge->flow;     /* FLOP */
}, GRB_GETRES_STR)

JIT_STR(void MF_extractMatrixFlow(double* flow, const MF_flowEdge* edge){*flow = edge->flow;}, GRB_EMFLOW_STR)

  
#endif
  
//------------------------------------------------------------------------------
// LAGraph_MaxFlow
//------------------------------------------------------------------------------

int LAGr_MaxFlow(double* f, GrB_Matrix* flow_mtx, LAGraph_Graph G, GrB_Index src, GrB_Index sink, char *msg){

#if LG_SUITESPARSE_GRAPHBLAS_V10

  //types
  GrB_Type GrB_FlowEdge = NULL ;
  GrB_Type GrB_ResultTuple = NULL ;
  GrB_Type GrB_CompareTuple = NULL ;

  GrB_Vector lvl = NULL ;
  GrB_UnaryOp GrB_GetResidual = NULL ;
  GrB_Matrix res_mat = NULL, res_matT = NULL ;
  LAGraph_Graph res_graph = NULL ;
    

  //to create R
  GrB_UnaryOp GrB_CreateResidualForward = NULL , GrB_CreateResidualBackward = NULL ;
  GrB_Matrix R = NULL ;

  //to init R with initial saturated flows
  GrB_Vector e = NULL, Re = NULL ;
  GrB_UnaryOp GrB_MakeFlow = NULL ;
  GrB_BinaryOp GrB_InitForwardFlows = NULL, GrB_InitBackwardFlows = NULL ;

  //create height vector
  GrB_Vector d = NULL ;

  //src and sink mask vec and n_active
  GrB_Vector src_and_sink = NULL ;
  GrB_Index n_active = INT64_MAX ;

  //semiring and vectors for y<e, struct> = R x d
  GrB_Vector y = NULL ;
  GrB_IndexUnaryOp GrB_Prune = NULL ;
  GxB_IndexBinaryOp GrB_RxdIndexMult = NULL ;
  GrB_BinaryOp GrB_RxdAdd = NULL, GrB_RxdMult = NULL ;
  GrB_Monoid GrB_RxdAddMonoid = NULL ;
  GrB_Semiring GrB_RxdSemiring = NULL ;
  GrB_Scalar theta = NULL ;
 
  //binary op and yd
  GrB_Vector yd = NULL ;
  GrB_BinaryOp GrB_CreateCompareVec = NULL ;
  
  //utility vectors, Matrix, and ops for mapping
  GrB_Matrix map = NULL ;
  GrB_Vector Jvec = NULL ;
  GrB_UnaryOp GrB_extractJ = NULL, GrB_extractYJ = NULL ;
  
  //map x e semiring
  GrB_Semiring GrB_MxeSemiring = NULL ;
  GrB_Monoid GrB_MxeAddMonoid = NULL ;
  GrB_BinaryOp GrB_MxeAdd = NULL, GrB_MxeMult = NULL ;
  GxB_IndexBinaryOp GrB_MxeIndexMult = NULL ;

  //residual flow vec
  GrB_Vector residual_vec = NULL ;
  GrB_UnaryOp GrB_extractFlows = NULL ;
  GrB_UnaryOp GrB_extractMatrixFlows = NULL ;
 
  //delta structures
  GrB_Vector delta_vec = NULL ;
  GrB_Matrix delta = NULL , delta_mat = NULL ;
 
  //update height
  GrB_BinaryOp GrB_UpdateHeight = NULL ;

  //update R structure
  GrB_BinaryOp GrB_UpdateFlows = NULL ;

  //scalars
  GrB_Scalar zero = NULL ;
  GrB_Scalar empty = NULL ;

  //invariant
  GrB_Vector invariant = NULL ;
  GrB_BinaryOp GrB_InvariantCheck = NULL ;
  GrB_Scalar check = NULL ;
  bool check_raw;

  //descriptor and matrix building
  GrB_Descriptor extract_desc = NULL ;

  //do input checks
  LG_TRY(LAGraph_CheckGraph(G, msg));
  LG_ASSERT (f != NULL, GrB_NULL_POINTER) ;
  (*f) = 0;
  GrB_Index nrows, n;
  GRB_TRY(GrB_Matrix_ncols(&n, G->A));
  GRB_TRY(GrB_Matrix_nrows(&nrows, G->A));
  LG_ASSERT_MSG(nrows == n, GrB_INVALID_VALUE, "Matrix must be square"); 
  LG_ASSERT_MSG(src < n && src >= 0 && sink < n && sink >= 0,
		GrB_INVALID_VALUE, "src and sink must be a value between [0, n)");
  LG_ASSERT_MSG(G->emin > 0, GrB_INVALID_VALUE,
		"the edge weights (capacities) must be greater than 0");

  //get adjacency matrix and its transpose
  GrB_Matrix A = G->A; 
  GrB_Matrix AT = NULL ;
  if (G->kind == LAGraph_ADJACENCY_UNDIRECTED)
  {
    // G is undirected, so A and AT are the same
    AT = G->A ;
  }
  else
  {
    // G is directed; get G->AT, which must be present
    AT = G->AT ;
    LG_ASSERT_MSG (AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required") ;
  }

  //create types for computation
  GRB_TRY(GxB_Type_new(&GrB_FlowEdge, sizeof(MF_flowEdge),
			 "MF_flowEdge", GRB_FLOWEDGE_STR));

  GRB_TRY(GxB_UnaryOp_new(&GrB_GetResidual, F_UNARY(MF_getResidual),
			    GrB_FP64, GrB_FlowEdge, "MF_getResidual",
			  GRB_GETRES_STR));
  GRB_TRY(GrB_Scalar_new(&check, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement(check, false));

  //create R
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidualForward,
			  F_UNARY(MF_CreateResidualForward),
			  GrB_FlowEdge , GrB_FP64,
			  "MF_CreateResidualForward", GRB_CRF_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidualBackward,
			  F_UNARY(MF_CreateResidualBackward),
			  GrB_FlowEdge , GrB_FP64,
			  "MF_CreateResidualBackward", GRB_CRB_STR));
  GRB_TRY(GrB_Matrix_new(&R, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_apply(R, NULL, NULL, GrB_CreateResidualForward, A, NULL));
  GRB_TRY(GrB_apply(R, A, NULL, GrB_CreateResidualBackward, AT, GrB_DESC_SC));

  //init R with initial saturated flows
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitForwardFlows,
			   F_BINARY(MF_initForwardFlows),
			   GrB_FlowEdge, GrB_FlowEdge,
			   GrB_FlowEdge, "MF_initForwardFlows", GRB_INITFLOWF_STR));
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitBackwardFlows,
			   F_BINARY(MF_initBackwardFlows),
			   GrB_FlowEdge, GrB_FlowEdge, GrB_FlowEdge,
			   "MF_initBackwardFlows", GRB_INITFLOWB_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_MakeFlow, F_UNARY(MF_MakeFlow),
			  GrB_FlowEdge, GrB_FP64, "MF_MakeFlow", GRB_MAKEF_STR));
  GRB_TRY(GrB_Vector_new(&Re, GrB_FlowEdge, n));
  GRB_TRY(GrB_Vector_new(&e, GrB_FP64, n));

  //extract n_active from e masking sink and src then assign to e
  GRB_TRY(GrB_Vector_new(&src_and_sink, GrB_BOOL, n));
  GRB_TRY (GrB_Vector_setElement (src_and_sink, true, sink)) ;
  GRB_TRY (GrB_Vector_setElement (src_and_sink, true, src)) ;

  
  //create flow vec
  GRB_TRY(GrB_Vector_new(&residual_vec, GrB_FP64, n));

  GRB_TRY(GrB_Matrix_new(&delta_mat, GrB_FP64, n, n));
  GRB_TRY(GrB_Matrix_new(&delta, GrB_FP64, n, n));
  GRB_TRY(GrB_Vector_new(&delta_vec, GrB_FP64, n));

  //update R structure
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateFlows,
			   F_BINARY(MF_updateFlow),
			   GrB_FlowEdge, GrB_FlowEdge,
			   GrB_FP64, "MF_updateFlow",
			   GRB_UPDATEFLOWS_STR));

  //create scalars
  GRB_TRY(GrB_Scalar_new(&zero, GrB_FP64));
  GRB_TRY(GrB_Scalar_setElement(zero, 0));
  GRB_TRY(GrB_Scalar_new (&empty, GrB_FP64)) ;

  GRB_TRY(GxB_UnaryOp_new(&GrB_extractMatrixFlows,
			  F_UNARY(MF_extractMatrixFlow),
			  GrB_FP64, GrB_FlowEdge, "MF_extractMatrixFlow", GRB_EMFLOW_STR));

  #ifdef COVERAGE
  // Just for test coverage, use 64-bit ints for n > 100.  Do not use this
  // rule in production!
  #define NBIG 100
  #else
  // For production use: 64-bit integers if n > 2^31
  #define NBIG INT32_MAX
  #endif
  if (n > NBIG){
    //create types for computation
    GRB_TRY(GxB_Type_new(&GrB_ResultTuple, sizeof(MF_resultTuple64),
			 "MF_resultTuple64", GRB_RESULTTUPLE_STR64));
    GRB_TRY(GxB_Type_new(&GrB_CompareTuple, sizeof(MF_compareTuple64),
			 "MF_compareTuple64", GRB_COMPARETUPLE_STR64));

    //invariant check
    #ifdef DBG
    GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
    GRB_TRY(GxB_BinaryOp_new(&GrB_InvariantCheck, F_BINARY(MF_CheckInvariant64),
			     GrB_BOOL, GrB_INT64,
			     GrB_ResultTuple, "MF_CheckInvariant64", GRB_INV_STR64));
    #endif

    //create and init d vector
    GRB_TRY(GrB_Vector_new(&d, GrB_INT64, n));
    GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(d, NULL, NULL, n, &src, 1, NULL));

    GRB_TRY(GxB_UnaryOp_new(&GrB_extractFlows,
			  F_UNARY(MF_extractFlow64),
			  GrB_FP64, GrB_ResultTuple,
			  "MF_extractFlow64", GRB_EXTRACTFLOW_STR64));

    //create semiring and vectors for y<e, struct> = R x d
    GRB_TRY(GrB_Scalar_new(&theta, GrB_INT64));
    GRB_TRY(GrB_Scalar_setElement_INT64(theta, 0));
    GRB_TRY(GrB_Vector_new(&y, GrB_ResultTuple, n));
    GRB_TRY(GxB_IndexBinaryOp_new(&GrB_RxdIndexMult,
				  F_INDEX_BINARY(MF_RxdMult64),
				  GrB_ResultTuple, GrB_FlowEdge,
				  GrB_INT64, GrB_INT64, "MF_RxdMult64",
				  GRB_RXDMULT_STR64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&GrB_RxdMult, GrB_RxdIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&GrB_RxdAdd, F_BINARY(MF_RxdAdd64),
			     GrB_ResultTuple, GrB_ResultTuple,
			     GrB_ResultTuple, "MF_RxdAdd64", GRB_RXDADD_STR64));
    MF_resultTuple64 id = {.d = INT64_MAX, .j = -1, .residual = 0};
    GRB_TRY(GrB_Monoid_new_UDT(&GrB_RxdAddMonoid, GrB_RxdAdd, &id));
    
    //create binary op and yd
    GRB_TRY(GrB_Vector_new(&yd, GrB_CompareTuple, n));
    GRB_TRY(GxB_BinaryOp_new(&GrB_CreateCompareVec, F_BINARY(MF_CreateCompareVec64),
			     GrB_CompareTuple, GrB_ResultTuple,
			     GrB_INT64, "MF_CreateCompareVec64",
			     GRB_CREATECOMPVEC_STR64));
    GRB_TRY(GxB_IndexUnaryOp_new(&GrB_Prune, (GxB_index_unary_function)
				 MF_Prune64, GrB_BOOL, GrB_ResultTuple,
				 GrB_INT64, "MF_Prune64", GRB_PRUNE_STR64));

    //create utility vectors, Matrix, and ops for mapping
    GRB_TRY(GrB_Vector_new(&Jvec, GrB_INT64, n));
    GRB_TRY(GxB_UnaryOp_new(&GrB_extractJ,
			    F_UNARY(MF_extractJ64),
			    GrB_INT64, GrB_CompareTuple,
			    "MF_extractJ64", GRB_EXTRACTJ_STR64));
    GRB_TRY(GxB_UnaryOp_new(&GrB_extractYJ, F_UNARY(MF_extractYJ64),
			    GrB_INT64, GrB_ResultTuple, "MF_extractYJ64",
			    GRB_EXTRACTYJ_STR64));

    //create map x e semiring
    GRB_TRY(GxB_IndexBinaryOp_new(&GrB_MxeIndexMult,
				  F_INDEX_BINARY(MF_MxeMult64),
				  GrB_ResultTuple, GrB_CompareTuple,
				  GrB_FP64, GrB_INT64, "MF_MxeMult64",
				  GRB_MXEMULT_STR64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&GrB_MxeMult, GrB_MxeIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&GrB_MxeAdd, F_BINARY(MF_MxeAdd64),
			     GrB_ResultTuple, GrB_ResultTuple,
			     GrB_ResultTuple, "MF_MxeAdd64",
			     GRB_MXEADD_STR64));
    GRB_TRY(GrB_Monoid_new_UDT(&GrB_MxeAddMonoid, GrB_MxeAdd, &id));

    //update height binary op
    GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateHeight,
			     F_BINARY(MF_updateHeight64),
			     GrB_INT64, GrB_INT64,
			     GrB_ResultTuple,
			     "MF_updateHeight64", GRB_UPDATEHEIGHT_STR64));
  }else{
    //create types for computation
    GRB_TRY(GxB_Type_new(&GrB_ResultTuple, sizeof(MF_resultTuple32),
			 "MF_resultTuple32", GRB_RESULTTUPLE_STR32));
    GRB_TRY(GxB_Type_new(&GrB_CompareTuple, sizeof(MF_compareTuple32),
			 "MF_compareTuple32", GRB_COMPARETUPLE_STR32));

    //invariant check
    #ifdef DBG
    GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
    GRB_TRY(GxB_BinaryOp_new(&GrB_InvariantCheck, F_BINARY(MF_CheckInvariant32),
			     GrB_BOOL, GrB_INT32,
			     GrB_ResultTuple, "MF_CheckInvariant32", GRB_INV_STR32));
    #endif

    GRB_TRY(GxB_UnaryOp_new(&GrB_extractFlows,
			  F_UNARY(MF_extractFlow32),
			  GrB_FP64, GrB_ResultTuple,
			  "MF_extractFlow32", GRB_EXTRACTFLOW_STR32));

    //create and init d vector
    GRB_TRY(GrB_Vector_new(&d, GrB_INT32, n));
    GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(d, NULL, NULL, n, &src, 1, NULL));
    
    
    //create semiring and vectors for y<e, struct> = R x d
    GRB_TRY(GrB_Scalar_new(&theta, GrB_INT32));
    GRB_TRY(GrB_Scalar_setElement_INT32(theta, 0));
    GRB_TRY(GrB_Vector_new(&y, GrB_ResultTuple, n));
    GRB_TRY(GxB_IndexBinaryOp_new(&GrB_RxdIndexMult,
				  F_INDEX_BINARY(MF_RxdMult32),
				  GrB_ResultTuple, GrB_FlowEdge,
				  GrB_INT32, GrB_INT32, "MF_RxdMult32",
				  GRB_RXDMULT_STR32));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&GrB_RxdMult, GrB_RxdIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&GrB_RxdAdd, F_BINARY(MF_RxdAdd32),
			     GrB_ResultTuple, GrB_ResultTuple,
			     GrB_ResultTuple, "MF_RxdAdd32", GRB_RXDADD_STR32));
    MF_resultTuple32 id = {.d = INT32_MAX, .j = -1, .residual = 0};

    GRB_TRY(GrB_Monoid_new_UDT(&GrB_RxdAddMonoid, GrB_RxdAdd, &id));

    //create binary op and yd
    GRB_TRY(GrB_Vector_new(&yd, GrB_CompareTuple, n));
    GRB_TRY(GxB_BinaryOp_new(&GrB_CreateCompareVec, F_BINARY(MF_CreateCompareVec32),
			     GrB_CompareTuple, GrB_ResultTuple,
			     GrB_INT32, "MF_CreateCompareVec32",
			     GRB_CREATECOMPVEC_STR32));
    GRB_TRY(GxB_IndexUnaryOp_new(&GrB_Prune, (GxB_index_unary_function)
				 MF_Prune32, GrB_BOOL, GrB_ResultTuple,
				 GrB_INT32, "MF_Prune32", GRB_PRUNE_STR32));

    //create utility vectors, Matrix, and ops for mapping
    GRB_TRY(GrB_Vector_new(&Jvec, GrB_INT32, n));
    GRB_TRY(GxB_UnaryOp_new(&GrB_extractJ,
			    F_UNARY(MF_extractJ32),
			    GrB_INT32, GrB_CompareTuple,
			    "MF_extractJ32", GRB_EXTRACTJ_STR32));
    GRB_TRY(GxB_UnaryOp_new(&GrB_extractYJ, F_UNARY(MF_extractYJ32),
			    GrB_INT32, GrB_ResultTuple,
			    "MF_extractYJ32", GRB_EXTRACTYJ_STR32));

    //create map x e semiring
    GRB_TRY(GxB_IndexBinaryOp_new(&GrB_MxeIndexMult,
				  F_INDEX_BINARY(MF_MxeMult32),
				  GrB_ResultTuple, GrB_CompareTuple,
				  GrB_FP64, GrB_INT32, "MF_MxeMult32",
				  GRB_MXEMULT_STR32));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&GrB_MxeMult, GrB_MxeIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&GrB_MxeAdd, F_BINARY(MF_MxeAdd32),
			     GrB_ResultTuple, GrB_ResultTuple,
			     GrB_ResultTuple, "MF_MxeAdd32",
			     GRB_MXEADD_STR32));
    GRB_TRY(GrB_Monoid_new_UDT(&GrB_MxeAddMonoid, GrB_MxeAdd, &id));

   
    //update height binary op
    GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateHeight,
			     F_BINARY(MF_updateHeight32),
			     GrB_INT32, GrB_INT32,
			     GrB_ResultTuple,
			     "MF_updateHeight32", GRB_UPDATEHEIGHT_STR32));

  }
  GRB_TRY(GrB_Matrix_new(&map, GrB_CompareTuple, n,n));

  GRB_TRY(GrB_Semiring_new(&GrB_RxdSemiring, GrB_RxdAddMonoid, GrB_RxdMult));
  GRB_TRY(GrB_Semiring_new(&GrB_MxeSemiring, GrB_MxeAddMonoid, GrB_MxeMult));
  
  int64_t iter = 0;

  //Create extract arrays
  GRB_TRY(GrB_Descriptor_new(&extract_desc));
  GRB_TRY(GrB_set(extract_desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)); 
  
  while(n_active > 0){

    // global relabeling, for the first iteration, and every 12 iterations
    // after that
    if((iter % 12 == 0) && (flow_mtx == NULL || iter == 0))
    {
      GRB_TRY(GrB_Matrix_new(&res_matT, GrB_FP64, n, n));
      GRB_TRY(GrB_Matrix_new(&res_mat, GrB_FP64, n, n));
      GRB_TRY(GrB_apply(res_mat, NULL,
			NULL, GrB_GetResidual, R, NULL)) ;  
      GRB_TRY(GrB_select(res_mat, NULL,
			 NULL, GrB_VALUEGT_FP64, res_mat, 0, NULL)) ;  
      GRB_TRY(GrB_transpose(res_matT, NULL, NULL, res_mat, NULL));
      LG_TRY(LAGraph_New(&res_graph,
			 &res_matT, LAGraph_ADJACENCY_DIRECTED, msg));
      res_graph->AT = res_mat;
      res_mat = NULL ;
      LG_TRY(LAGraph_Cached_OutDegree(res_graph, msg));
      LG_TRY(LAGr_BreadthFirstSearch(&lvl,
				     NULL,
				     res_graph, sink, msg));

      // d<![src,sink],struct> = lvl
      GRB_TRY(GrB_assign(d, src_and_sink,
			 NULL, lvl, GrB_ALL, n, GrB_DESC_SC));
      // d<!lvl,struct> = n
      GRB_TRY(GrB_assign(d, lvl, NULL, n,
			 GrB_ALL, n, GrB_DESC_SC));

      if(iter == 0){
	GRB_TRY(GrB_extract(e, lvl, NULL, A, GrB_ALL, n, src, GrB_DESC_ST0));
	GRB_TRY(GrB_apply(Re, NULL, NULL, GrB_MakeFlow, e, NULL));
        GRB_TRY(GrB_assign(R, NULL, GrB_InitForwardFlows, Re, src, GrB_ALL, n, NULL));
        GRB_TRY(GrB_assign(R, NULL, GrB_InitBackwardFlows, Re, GrB_ALL, n, src, NULL));
	LG_TRY (LG_augment_maxflow (f, e, sink, src_and_sink, &n_active, msg)) ;
      }
      else{
	GrB_assign (e, lvl, NULL, empty, GrB_ALL, n, GrB_DESC_SC) ;
      }

      GrB_free(&lvl);
      LG_TRY(LAGraph_Delete(&res_graph, msg));
      GRB_TRY(GrB_Vector_nvals(&n_active, e));
      if(n_active == 0){
	break;
      }
    }
    
    GRB_TRY(GrB_mxv(y, e, NULL,
		    GrB_RxdSemiring, R, d, GrB_DESC_RS));
    GRB_TRY(GrB_select(y, NULL, NULL, GrB_Prune, y, -1, NULL));

    //create yd vector of type compare tuple
    GRB_TRY(GrB_eWiseMult(yd, NULL, NULL, GrB_CreateCompareVec, y,  d, NULL));

    //create map matrix from yd
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractJ, yd, NULL));
    GRB_TRY(GrB_Matrix_clear(map));
    GRB_TRY(GrB_Matrix_build(map, yd,
			     Jvec, yd, GxB_IGNORE_DUP, extract_desc));
    
    //make e dense for map computation
    GRB_TRY(GrB_assign(e, e, NULL, 0, GrB_ALL, n, GrB_DESC_SC));

    //y = map x e
    GRB_TRY(GrB_mxv(y, NULL, NULL, GrB_MxeSemiring, map, e, NULL));
    GRB_TRY(GrB_select(y, NULL, NULL, GrB_Prune, y, -1, NULL));
   
    //relable, update heights
    GRB_TRY(GrB_eWiseMult(d, y, NULL, GrB_UpdateHeight, d, y, GrB_DESC_S));

    #ifdef DBG
    //assert correct labels
        GRB_TRY(GrB_eWiseMult(invariant, y,
			      NULL, GrB_InvariantCheck, d, y, GrB_DESC_RS));
	GRB_TRY(GrB_reduce(check, NULL,
			   GrB_LAND_MONOID_BOOL, invariant, NULL));
	GRB_TRY(GrB_Scalar_extractElement(&check_raw, check));
	ASSERT(check_raw == true);
    #endif

    //extract residual flows from y
    GRB_TRY(GrB_apply(residual_vec, NULL, NULL, GrB_extractFlows, y, NULL));

    //.min(flow_vec and e)
    GRB_TRY(GrB_eWiseMult(delta_vec, NULL, NULL,
			  GrB_MIN_FP64, residual_vec, e, NULL));    /* FLOP */
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractYJ, y, NULL));
    GRB_TRY(GrB_Matrix_clear(delta));
    GRB_TRY(GxB_Matrix_build_Vector(delta, delta_vec,
				    Jvec, delta_vec, GxB_IGNORE_DUP, extract_desc));

    //make delta anti-symmetric: delta_mat = (delta - delta')
    GRB_TRY(GxB_eWiseUnion(delta_mat, NULL, NULL, GrB_MINUS_FP64,   /* FLOP */
			   delta, zero, delta, zero, GrB_DESC_T1));

    //update R
    // R<delta_mat> = UpdateFlows (R, delta_mat) using eWiseMult
    GRB_TRY(GrB_eWiseMult(R, delta_mat, NULL, GrB_UpdateFlows,  /* FLOP */
			  R, delta_mat, GrB_DESC_S));

    //reduce delta_mat to delta_vec
    // delta_vec = sum (delta_mat), summing up each row of delta_mat
    GRB_TRY(GrB_reduce(delta_vec, NULL, NULL,
		       GrB_PLUS_FP64, delta_mat, GrB_DESC_T0)); /* FLOP */

    //add to e
    // e<delta_vec> += delta_vec
    GRB_TRY(GrB_assign(e, delta_vec, GrB_PLUS_FP64, /* FLOP */
		       delta_vec, GrB_ALL, n, GrB_DESC_S));
    
    // augment maxflow for all active nodes
    LG_TRY (LG_augment_maxflow (f, e, sink, src_and_sink,
				&n_active, msg)) ;

    ++iter;
    
  }

  if (flow_mtx != NULL)
  {
    GRB_TRY(GrB_apply(*flow_mtx, NULL, NULL, GrB_extractMatrixFlows, R, NULL));
    GRB_TRY(GrB_select(*flow_mtx, NULL, NULL, GrB_VALUEGE_FP64, *flow_mtx, 0, NULL));
  }

  #ifdef COVERAGE
  // The GrB_MxeAdd operator is not tested via the call to GrB_mxv with the
  // GrB_MxeSemiring above, so test it via the GrB_MxeAddMonoid.
  GrB_free(&y);
  GRB_TRY(GrB_Vector_new(&y, GrB_ResultTuple, 3));
  if (n > NBIG)
  {
    MF_resultTuple64 a = {.d = 1, .j = 2, .residual = 3};
    MF_resultTuple64 b = {.d = 4, .j = 5, .residual = 6};
    GRB_TRY (GrB_Vector_setElement_UDT (y, (void *) &a, 0)) ;
    GRB_TRY (GrB_Vector_setElement_UDT (y, (void *) &b, 0)) ;
    MF_resultTuple64 c = {.d = 0, .j = 0, .residual = 0};
    GRB_TRY (GrB_Vector_reduce_UDT ((void *) &c, NULL, GrB_MxeAddMonoid, y, NULL)) ;
//  printf ("c: resid %g j %g d %g\n", c.residual, (double) c.j, (double) c.d) ;
    LG_ASSERT ((c.residual == 6 && c.j == 5 && c.d == 4), GrB_PANIC) ;
  }
  else
  {
    MF_resultTuple32 a = {.d = 1, .j = 2, .residual = 3};
    MF_resultTuple32 b = {.d = 4, .j = 5, .residual = 6};
    GRB_TRY (GrB_Vector_setElement_UDT (y, (void *) &a, 0)) ;
    GRB_TRY (GrB_Vector_setElement_UDT (y, (void *) &b, 0)) ;
    MF_resultTuple32 c = {.d = 0, .j = 0, .residual = 0};
    GRB_TRY (GrB_Vector_reduce_UDT ((void *) &c, NULL, GrB_MxeAddMonoid, y, NULL)) ;
//  printf ("c: resid %g j %g d %g\n", c.residual, (double) c.j, (double) c.d) ;
    LG_ASSERT ((c.residual == 6 && c.j == 5 && c.d == 4), GrB_PANIC) ;
  }
  #endif

  LG_FREE_ALL;
  return GrB_SUCCESS;
#else
  return GrB_NOT_IMPLEMENTED ;
#endif
}
