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
    GrB_free(&FlowEdge);                                                       \
    GrB_free(&CompareTuple);                                                   \
    GrB_free(&ResultTuple);                                                    \
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
    GrB_free(&Prune);                                                          \
    GrB_free(&UpdateFlows);                                                    \
    GrB_free(&UpdateHeight);                                                   \
    GrB_free(&extractFlows);                                                   \
    GrB_free(&MxeIndexMult);                                                   \
    GrB_free(&MxeMult);                                                        \
    GrB_free(&MxeAdd);                                                         \
    GrB_free(&MxeAddMonoid);                                                   \
    GrB_free(&MxeSemiring);                                                    \
    GrB_free(&extractJ);                                                       \
    GrB_free(&CreateCompareVec);                                               \
    GrB_free(&RxdSemiring);                                                    \
    GrB_free(&RxdAdd);                                                         \
    GrB_free(&RxdAddMonoid);                                                   \
    GrB_free(&RxdIndexMult);                                                   \
    GrB_free(&RxdMult);                                                        \
    GrB_free(&InitForwardFlows);                                               \
    GrB_free(&InitBackwardFlows);                                              \
    GrB_free(&CreateResidualForward);                                          \
    GrB_free(&CreateResidualBackward);                                         \
    GrB_free(&zero);                                                           \
    GrB_free(&empty);                                                          \
    GrB_free(&Re);                                                             \
    GrB_free(&invariant);                                                      \
    GrB_free(&InvariantCheck);                                                 \
    GrB_free(&check);                                                          \
    GrB_free(&extractYJ);                                                      \
    GrB_free(&extract_desc);                                                   \
    GrB_free(&residual_vec);						       \
    GrB_free(&MakeFlow);						       \
    GrB_free(&GetResidual);				                       \
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
  } MF_flowEdge;, FLOWEDGE_STR)

JIT_STR(typedef struct{
  double residual;
  int64_t j;
  int64_t d;
  } MF_resultTuple64;, RESULTTUPLE_STR64)

JIT_STR(typedef struct{
  double residual;
  int32_t j;
  int32_t d;
  } MF_resultTuple32;, RESULTTUPLE_STR32)


JIT_STR(typedef struct{
  double residual;
  int64_t di;
  int64_t y_dmin;
  int64_t j;
  } MF_compareTuple64;, COMPARETUPLE_STR64)

JIT_STR(typedef struct{
  double residual;
  int32_t di;
  int32_t y_dmin;
  int32_t j;
  int32_t unused;   /* to pad the struct to 24 bytes */
  } MF_compareTuple32;, COMPARETUPLE_STR32) // 24 bytes: padded


JIT_STR(void MF_CreateResidualForward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = (*y);
  }, CRF_STR)

JIT_STR(void MF_CreateResidualBackward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = 0;
  }, CRB_STR)


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
  }, RXDMULT_STR64)


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
}, RXDMULT_STR32)

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
  }, RXDADD_STR64)

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
  }, RXDADD_STR32)



JIT_STR(void MF_extractFlow64(double *z, const MF_resultTuple64 *x)
	{ (*z) = x->residual; }, EXTRACTFLOW_STR64)

JIT_STR(void MF_extractFlow32(double *z, const MF_resultTuple32 *x)
	{ (*z) = x->residual; }, EXTRACTFLOW_STR32)


JIT_STR(void MF_updateFlow(MF_flowEdge *z,
			   const MF_flowEdge *x, const double *y) {
  z->capacity = x->capacity;
  z->flow = x->flow + (*y);
  }, UPDATEFLOWS_STR)


JIT_STR(void MF_updateHeight64(int64_t *z,
			     const int64_t *x, const MF_resultTuple64 *y) {
  if((*x) < y->d+1){
    (*z) = y->d + 1;
  }
  else {
    (*z) = (*x);
  }
  }, UPDATEHEIGHT_STR64)

JIT_STR(void MF_updateHeight32(int32_t *z,
			     const int32_t *x, const MF_resultTuple32 *y) {
  if((*x) < y->d+1){
    (*z) = y->d + 1;
  }
  else {
    (*z) = (*x);
  }
  }, UPDATEHEIGHT_STR32)


JIT_STR(void MF_extractJ64(int64_t *z, const MF_compareTuple64 *x) { (*z) = x->j; }, EXTRACTJ_STR64)

JIT_STR(void MF_extractJ32(int32_t *z, const MF_compareTuple32 *x) { (*z) = x->j; }, EXTRACTJ_STR32)


JIT_STR(void MF_extractYJ64(int64_t *z, const MF_resultTuple64 *x) {
  (*z) = x->j;
  }, EXTRACTYJ_STR64)

JIT_STR(void MF_extractYJ32(int32_t *z, const MF_resultTuple32 *x) {
  (*z) = x->j;
  }, EXTRACTYJ_STR32)


JIT_STR(void MF_initForwardFlows(MF_flowEdge * z,
				 const MF_flowEdge * x, const MF_flowEdge * y){
  z->flow = y->flow + x->flow;
  z->capacity = x->capacity;
  }, INITFLOWF_STR)


JIT_STR(void MF_initBackwardFlows(MF_flowEdge * z,
				  const MF_flowEdge * x, const MF_flowEdge * y){
  z->flow = x->flow - y->flow;
  z->capacity = x->capacity;
  }, INITFLOWB_STR)

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
  }, MXEMULT_STR64)

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
  }, MXEMULT_STR32)


// Note: the additive monoid is not actually used in the call to GrB_mxv below,
// because any given node only pushes to one neighbor at a time.  As a result,
// no reduction is needed in GrB_mxv.  The semiring still needs a monoid,
// however.
JIT_STR(void MF_MxeAdd64(MF_resultTuple64 * z,
		       const MF_resultTuple64 * x, const MF_resultTuple64 * y){
    (*z) = (*y) ;
  }, MXEADD_STR64)

JIT_STR(void MF_MxeAdd32(MF_resultTuple32 * z,
		       const MF_resultTuple32 * x, const MF_resultTuple32 * y){
    (*z) = (*y) ;
  }, MXEADD_STR32)


JIT_STR(void MF_CreateCompareVec64(MF_compareTuple64 *comp,
				 const MF_resultTuple64 *res, const int64_t *height) {
  comp->di = (*height);
  comp->j = res->j;
  comp->residual = res->residual;
  comp->y_dmin = res->d;
  }, CREATECOMPVEC_STR64)

JIT_STR(void MF_CreateCompareVec32(MF_compareTuple32 *comp,
				 const MF_resultTuple32 *res, const int32_t *height) {
  comp->di = (*height);
  comp->j = res->j;
  comp->residual = res->residual;
  comp->y_dmin = res->d;
  comp->unused = 0 ;
  }, CREATECOMPVEC_STR32)


JIT_STR(void MF_Prune64(bool * z, const MF_resultTuple64 * x,
		      GrB_Index ix, GrB_Index jx, const int64_t * theta){
  *z = (x->j != *theta) ;
  }, PRUNE_STR64)

JIT_STR(void MF_Prune32(bool * z, const MF_resultTuple32 * x,
		      GrB_Index ix, GrB_Index jx, const int32_t * theta){
  *z = (x->j != *theta) ;
  }, PRUNE_STR32)


JIT_STR(void MF_MakeFlow(MF_flowEdge * flow_edge, const double * flow){
  flow_edge->capacity = 0;
  flow_edge->flow = (*flow);
  }, MAKEF_STR)

#ifdef DBG
JIT_STR(void MF_CheckInvariant64(bool *z, const int64_t *height,
			       const MF_resultTuple64 *result) {
  (*z) = ((*height) == result->d+1);
  }, INV_STR64)

JIT_STR(void MF_CheckInvariant32(bool *z, const int32_t *height,
			       const MF_resultTuple32 *result) {
  (*z) = ((*height) == result->d+1);
  }, INV_STR32)
#endif

JIT_STR(void MF_getResidual(double * res, const MF_flowEdge * flow_edge){
    (*res) = flow_edge->capacity - flow_edge->flow;     /* FLOP */
}, GETRES_STR)

JIT_STR(void MF_extractMatrixFlow(double* flow, const MF_flowEdge* edge){*flow = edge->flow;}, EMFLOW_STR)

#endif

//------------------------------------------------------------------------------
// LAGraph_MaxFlow
//------------------------------------------------------------------------------

int LAGr_MaxFlow(double* f, GrB_Matrix* flow_mtx, LAGraph_Graph G, GrB_Index src, GrB_Index sink, char *msg){

#if LG_SUITESPARSE_GRAPHBLAS_V10

  //----------------------------------------------------------------------------
  // declare variables
  //----------------------------------------------------------------------------

  //types
  GrB_Type FlowEdge = NULL ;
  GrB_Type ResultTuple = NULL ;
  GrB_Type CompareTuple = NULL ;

  GrB_Vector lvl = NULL ;
  GrB_UnaryOp GetResidual = NULL ;
  GrB_Matrix res_mat = NULL, res_matT = NULL ;
  LAGraph_Graph res_graph = NULL ;

  //to create R
  GrB_UnaryOp CreateResidualForward = NULL, CreateResidualBackward = NULL ;
  GrB_Matrix R = NULL ;

  //to init R with initial saturated flows
  GrB_Vector e = NULL, Re = NULL ;
  GrB_UnaryOp MakeFlow = NULL ;
  GrB_BinaryOp InitForwardFlows = NULL, InitBackwardFlows = NULL ;

  //create height vector
  GrB_Vector d = NULL ;

  //src and sink mask vec and n_active
  GrB_Vector src_and_sink = NULL ;
  GrB_Index n_active = INT64_MAX ;

  //semiring and vectors for y<e, struct> = R x d
  GrB_Vector y = NULL ;
  GrB_IndexUnaryOp Prune = NULL ;
  GxB_IndexBinaryOp RxdIndexMult = NULL ;
  GrB_BinaryOp RxdAdd = NULL, RxdMult = NULL ;
  GrB_Monoid RxdAddMonoid = NULL ;
  GrB_Semiring RxdSemiring = NULL ;
  GrB_Scalar theta = NULL ;

  //binary op and yd
  GrB_Vector yd = NULL ;
  GrB_BinaryOp CreateCompareVec = NULL ;

  //utility vectors, Matrix, and ops for mapping
  GrB_Matrix map = NULL ;
  GrB_Vector Jvec = NULL ;
  GrB_UnaryOp extractJ = NULL, extractYJ = NULL ;

  //map x e semiring
  GrB_Semiring MxeSemiring = NULL ;
  GrB_Monoid MxeAddMonoid = NULL ;
  GrB_BinaryOp MxeAdd = NULL, MxeMult = NULL ;
  GxB_IndexBinaryOp MxeIndexMult = NULL ;

  //residual flow vec
  GrB_Vector residual_vec = NULL ;
  GrB_UnaryOp extractFlows = NULL ;
  GrB_UnaryOp extractMatrixFlows = NULL ;

  //delta structures
  GrB_Vector delta_vec = NULL ;
  GrB_Matrix delta = NULL , delta_mat = NULL ;

  //update height
  GrB_BinaryOp UpdateHeight = NULL ;

  //update R structure
  GrB_BinaryOp UpdateFlows = NULL ;

  //scalars
  GrB_Scalar zero = NULL ;
  GrB_Scalar empty = NULL ;

  //invariant
  GrB_Vector invariant = NULL ;
  GrB_BinaryOp InvariantCheck = NULL ;
  GrB_Scalar check = NULL ;
  bool check_raw;

  //descriptor and matrix building
  GrB_Descriptor extract_desc = NULL ;

  //----------------------------------------------------------------------------
  // check inputs
  //----------------------------------------------------------------------------

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

  //----------------------------------------------------------------------------
  // create types, operators, matrices, and vectors
  //----------------------------------------------------------------------------

  //create types for computation
  GRB_TRY(GxB_Type_new(&FlowEdge, sizeof(MF_flowEdge),
			 "MF_flowEdge", FLOWEDGE_STR));

  GRB_TRY(GxB_UnaryOp_new(&GetResidual, F_UNARY(MF_getResidual),
			  GrB_FP64, FlowEdge, "MF_getResidual",
			  GETRES_STR));
  GRB_TRY(GrB_Scalar_new(&check, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement(check, false));

  //create R
  GRB_TRY(GxB_UnaryOp_new(&CreateResidualForward,
			  F_UNARY(MF_CreateResidualForward),
			  FlowEdge , GrB_FP64,
			  "MF_CreateResidualForward", CRF_STR));
  GRB_TRY(GxB_UnaryOp_new(&CreateResidualBackward,
			  F_UNARY(MF_CreateResidualBackward),
			  FlowEdge , GrB_FP64,
			  "MF_CreateResidualBackward", CRB_STR));
  GRB_TRY(GrB_Matrix_new(&R, FlowEdge, n, n));
  GRB_TRY(GrB_apply(R, NULL, NULL, CreateResidualForward, A, NULL));
  GRB_TRY(GrB_apply(R, A, NULL, CreateResidualBackward, AT, GrB_DESC_SC));

  //init R with initial saturated flows
  GRB_TRY(GxB_BinaryOp_new(&InitForwardFlows,
			   F_BINARY(MF_initForwardFlows),
			   FlowEdge, FlowEdge, FlowEdge,
                           "MF_initForwardFlows", INITFLOWF_STR));
  GRB_TRY(GxB_BinaryOp_new(&InitBackwardFlows,
			   F_BINARY(MF_initBackwardFlows),
			   FlowEdge, FlowEdge, FlowEdge,
			   "MF_initBackwardFlows", INITFLOWB_STR));
  GRB_TRY(GxB_UnaryOp_new(&MakeFlow, F_UNARY(MF_MakeFlow),
			  FlowEdge, GrB_FP64, "MF_MakeFlow", MAKEF_STR));
  GRB_TRY(GrB_Vector_new(&Re, FlowEdge, n));
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
  GRB_TRY(GxB_BinaryOp_new(&UpdateFlows, F_BINARY(MF_updateFlow),
			   FlowEdge, FlowEdge, GrB_FP64, "MF_updateFlow",
			   UPDATEFLOWS_STR));

  //create scalars
  GRB_TRY(GrB_Scalar_new(&zero, GrB_FP64));
  GRB_TRY(GrB_Scalar_setElement(zero, 0));
  GRB_TRY(GrB_Scalar_new (&empty, GrB_FP64)) ;

  GRB_TRY(GxB_UnaryOp_new(&extractMatrixFlows,
			  F_UNARY(MF_extractMatrixFlow), GrB_FP64, FlowEdge,
                          "MF_extractMatrixFlow", EMFLOW_STR));

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
    GRB_TRY(GxB_Type_new(&ResultTuple, sizeof(MF_resultTuple64),
			 "MF_resultTuple64", RESULTTUPLE_STR64));
    GRB_TRY(GxB_Type_new(&CompareTuple, sizeof(MF_compareTuple64),
			 "MF_compareTuple64", COMPARETUPLE_STR64));

    //invariant check
    #ifdef DBG
    GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
    GRB_TRY(GxB_BinaryOp_new(&InvariantCheck, F_BINARY(MF_CheckInvariant64),
			     GrB_BOOL, GrB_INT64, ResultTuple,
                             "MF_CheckInvariant64", INV_STR64));
    #endif

    //create and init d vector
    GRB_TRY(GrB_Vector_new(&d, GrB_INT64, n));
    GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(d, NULL, NULL, n, &src, 1, NULL));

    GRB_TRY(GxB_UnaryOp_new(&extractFlows, F_UNARY(MF_extractFlow64),
			  GrB_FP64, ResultTuple,
			  "MF_extractFlow64", EXTRACTFLOW_STR64));

    //create semiring and vectors for y<e, struct> = R x d
    GRB_TRY(GrB_Scalar_new(&theta, GrB_INT64));
    GRB_TRY(GrB_Scalar_setElement_INT64(theta, 0));
    GRB_TRY(GrB_Vector_new(&y, ResultTuple, n));
    GRB_TRY(GxB_IndexBinaryOp_new(&RxdIndexMult, F_INDEX_BINARY(MF_RxdMult64),
				  ResultTuple, FlowEdge, GrB_INT64, GrB_INT64,
                                  "MF_RxdMult64", RXDMULT_STR64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&RxdMult, RxdIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&RxdAdd, F_BINARY(MF_RxdAdd64),
			     ResultTuple, ResultTuple, ResultTuple,
                             "MF_RxdAdd64", RXDADD_STR64));
    MF_resultTuple64 id = {.d = INT64_MAX, .j = -1, .residual = 0};
    GRB_TRY(GrB_Monoid_new_UDT(&RxdAddMonoid, RxdAdd, &id));

    //create binary op and yd
    GRB_TRY(GrB_Vector_new(&yd, CompareTuple, n));
    GRB_TRY(GxB_BinaryOp_new(&CreateCompareVec, F_BINARY(MF_CreateCompareVec64),
			     CompareTuple, ResultTuple, GrB_INT64,
                             "MF_CreateCompareVec64", CREATECOMPVEC_STR64));
    GRB_TRY(GxB_IndexUnaryOp_new(&Prune, (GxB_index_unary_function)
				 MF_Prune64, GrB_BOOL, ResultTuple,
				 GrB_INT64, "MF_Prune64", PRUNE_STR64));

    //create utility vectors, Matrix, and ops for mapping
    GRB_TRY(GrB_Vector_new(&Jvec, GrB_INT64, n));
    GRB_TRY(GxB_UnaryOp_new(&extractJ, F_UNARY(MF_extractJ64),
			    GrB_INT64, CompareTuple,
			    "MF_extractJ64", EXTRACTJ_STR64));
    GRB_TRY(GxB_UnaryOp_new(&extractYJ, F_UNARY(MF_extractYJ64),
			    GrB_INT64, ResultTuple, "MF_extractYJ64",
			    EXTRACTYJ_STR64));

    //create map x e semiring
    GRB_TRY(GxB_IndexBinaryOp_new(&MxeIndexMult, F_INDEX_BINARY(MF_MxeMult64),
				  ResultTuple, CompareTuple,
				  GrB_FP64, GrB_INT64, "MF_MxeMult64",
				  MXEMULT_STR64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MxeMult, MxeIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&MxeAdd, F_BINARY(MF_MxeAdd64),
			     ResultTuple, ResultTuple, ResultTuple,
                             "MF_MxeAdd64", MXEADD_STR64));
    GRB_TRY(GrB_Monoid_new_UDT(&MxeAddMonoid, MxeAdd, &id));

    //update height binary op
    GRB_TRY(GxB_BinaryOp_new(&UpdateHeight, F_BINARY(MF_updateHeight64),
			     GrB_INT64, GrB_INT64, ResultTuple,
			     "MF_updateHeight64", UPDATEHEIGHT_STR64));
  }else{
    //create types for computation
    GRB_TRY(GxB_Type_new(&ResultTuple, sizeof(MF_resultTuple32),
			 "MF_resultTuple32", RESULTTUPLE_STR32));
    GRB_TRY(GxB_Type_new(&CompareTuple, sizeof(MF_compareTuple32),
			 "MF_compareTuple32", COMPARETUPLE_STR32));

    //invariant check
    #ifdef DBG
    GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
    GRB_TRY(GxB_BinaryOp_new(&InvariantCheck, F_BINARY(MF_CheckInvariant32),
			     GrB_BOOL, GrB_INT32,
			     ResultTuple, "MF_CheckInvariant32", INV_STR32));
    #endif

    GRB_TRY(GxB_UnaryOp_new(&extractFlows, F_UNARY(MF_extractFlow32),
			  GrB_FP64, ResultTuple,
			  "MF_extractFlow32", EXTRACTFLOW_STR32));

    //create and init d vector
    GRB_TRY(GrB_Vector_new(&d, GrB_INT32, n));
    GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(d, NULL, NULL, n, &src, 1, NULL));

    //create semiring and vectors for y<e, struct> = R x d
    GRB_TRY(GrB_Scalar_new(&theta, GrB_INT32));
    GRB_TRY(GrB_Scalar_setElement_INT32(theta, 0));
    GRB_TRY(GrB_Vector_new(&y, ResultTuple, n));
    GRB_TRY(GxB_IndexBinaryOp_new(&RxdIndexMult, F_INDEX_BINARY(MF_RxdMult32),
				  ResultTuple, FlowEdge, GrB_INT32, GrB_INT32,
                                  "MF_RxdMult32", RXDMULT_STR32));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&RxdMult, RxdIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&RxdAdd, F_BINARY(MF_RxdAdd32),
			     ResultTuple, ResultTuple, ResultTuple,
                             "MF_RxdAdd32", RXDADD_STR32));
    MF_resultTuple32 id = {.d = INT32_MAX, .j = -1, .residual = 0};

    GRB_TRY(GrB_Monoid_new_UDT(&RxdAddMonoid, RxdAdd, &id));

    //create binary op and yd
    GRB_TRY(GrB_Vector_new(&yd, CompareTuple, n));
    GRB_TRY(GxB_BinaryOp_new(&CreateCompareVec, F_BINARY(MF_CreateCompareVec32),
			     CompareTuple, ResultTuple, GrB_INT32,
                             "MF_CreateCompareVec32", CREATECOMPVEC_STR32));
    GRB_TRY(GxB_IndexUnaryOp_new(&Prune, (GxB_index_unary_function)
				 MF_Prune32, GrB_BOOL, ResultTuple, GrB_INT32,
                                 "MF_Prune32", PRUNE_STR32));

    //create utility vectors, Matrix, and ops for mapping
    GRB_TRY(GrB_Vector_new(&Jvec, GrB_INT32, n));
    GRB_TRY(GxB_UnaryOp_new(&extractJ, F_UNARY(MF_extractJ32),
			    GrB_INT32, CompareTuple,
			    "MF_extractJ32", EXTRACTJ_STR32));
    GRB_TRY(GxB_UnaryOp_new(&extractYJ, F_UNARY(MF_extractYJ32),
			    GrB_INT32, ResultTuple,
			    "MF_extractYJ32", EXTRACTYJ_STR32));

    //create map x e semiring
    GRB_TRY(GxB_IndexBinaryOp_new(&MxeIndexMult, F_INDEX_BINARY(MF_MxeMult32),
				  ResultTuple, CompareTuple,
				  GrB_FP64, GrB_INT32, "MF_MxeMult32",
				  MXEMULT_STR32));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MxeMult, MxeIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&MxeAdd, F_BINARY(MF_MxeAdd32),
			     ResultTuple, ResultTuple, ResultTuple,
                             "MF_MxeAdd32", MXEADD_STR32));
    GRB_TRY(GrB_Monoid_new_UDT(&MxeAddMonoid, MxeAdd, &id));

    //update height binary op
    GRB_TRY(GxB_BinaryOp_new(&UpdateHeight, F_BINARY(MF_updateHeight32),
			     GrB_INT32, GrB_INT32, ResultTuple,
			     "MF_updateHeight32", UPDATEHEIGHT_STR32));

  }
  GRB_TRY(GrB_Matrix_new(&map, CompareTuple, n,n));

  GRB_TRY(GrB_Semiring_new(&RxdSemiring, RxdAddMonoid, RxdMult));
  GRB_TRY(GrB_Semiring_new(&MxeSemiring, MxeAddMonoid, MxeMult));

  int64_t iter = 0;

  //Create extract arrays
  GRB_TRY(GrB_Descriptor_new(&extract_desc));
  GRB_TRY(GrB_set(extract_desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST));

  //----------------------------------------------------------------------------
  // compute the max flow
  //----------------------------------------------------------------------------

  while(n_active > 0){

    // global relabeling, for the first iteration, and every 12 iterations
    // after that.  If the flow matrix is to be returned, global relabelling
    // can only be done on the first iteration.
    if((iter % 12 == 0) && (flow_mtx == NULL || iter == 0))
    {
      GRB_TRY(GrB_Matrix_new(&res_matT, GrB_FP64, n, n));
      GRB_TRY(GrB_Matrix_new(&res_mat, GrB_FP64, n, n));
      GRB_TRY(GrB_apply(res_mat, NULL, NULL, GetResidual, R, NULL)) ;
      GRB_TRY(GrB_select(res_mat, NULL, NULL, GrB_VALUEGT_FP64, res_mat, 0, NULL)) ;
      GRB_TRY(GrB_transpose(res_matT, NULL, NULL, res_mat, NULL));
      LG_TRY(LAGraph_New(&res_graph,
			 &res_matT, LAGraph_ADJACENCY_DIRECTED, msg));
      res_graph->AT = res_mat;
      res_mat = NULL ;
      LG_TRY(LAGraph_Cached_OutDegree(res_graph, msg));
      LG_TRY(LAGr_BreadthFirstSearch(&lvl, NULL, res_graph, sink, msg));

      // d<![src,sink],struct> = lvl
      GRB_TRY(GrB_assign(d, src_and_sink, NULL, lvl, GrB_ALL, n, GrB_DESC_SC));
      // d<!lvl,struct> = n
      GRB_TRY(GrB_assign(d, lvl, NULL, n, GrB_ALL, n, GrB_DESC_SC));

      if(iter == 0){
	GRB_TRY(GrB_extract(e, lvl, NULL, A, GrB_ALL, n, src, GrB_DESC_ST0));
	GRB_TRY(GrB_apply(Re, NULL, NULL, MakeFlow, e, NULL));
        GRB_TRY(GrB_assign(R, NULL, InitForwardFlows, Re, src, GrB_ALL, n, NULL));
        GRB_TRY(GrB_assign(R, NULL, InitBackwardFlows, Re, GrB_ALL, n, src, NULL));
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

    GRB_TRY(GrB_mxv(y, e, NULL, RxdSemiring, R, d, GrB_DESC_RS));
    GRB_TRY(GrB_select(y, NULL, NULL, Prune, y, -1, NULL));

    //create yd vector of type compare tuple
    GRB_TRY(GrB_eWiseMult(yd, NULL, NULL, CreateCompareVec, y,  d, NULL));

    //create map matrix from yd
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, extractJ, yd, NULL));
    GRB_TRY(GrB_Matrix_clear(map));
    GRB_TRY(GrB_Matrix_build(map, yd, Jvec, yd, GxB_IGNORE_DUP, extract_desc));

    //make e dense for map computation
    GRB_TRY(GrB_assign(e, e, NULL, 0, GrB_ALL, n, GrB_DESC_SC));

    //y = map x e
    GRB_TRY(GrB_mxv(y, NULL, NULL, MxeSemiring, map, e, NULL));
    GRB_TRY(GrB_select(y, NULL, NULL, Prune, y, -1, NULL));

    //relable, update heights
    GRB_TRY(GrB_eWiseMult(d, y, NULL, UpdateHeight, d, y, GrB_DESC_S));

    #ifdef DBG
    //assert correct labels
        GRB_TRY(GrB_eWiseMult(invariant, y,
			      NULL, InvariantCheck, d, y, GrB_DESC_RS));
	GRB_TRY(GrB_reduce(check, NULL,
			   GrB_LAND_MONOID_BOOL, invariant, NULL));
	GRB_TRY(GrB_Scalar_extractElement(&check_raw, check));
	ASSERT(check_raw == true);
    #endif

    //extract residual flows from y
    GRB_TRY(GrB_apply(residual_vec, NULL, NULL, extractFlows, y, NULL));

    //.min(flow_vec and e)
    GRB_TRY(GrB_eWiseMult(delta_vec, NULL, NULL,
			  GrB_MIN_FP64, residual_vec, e, NULL));    /* FLOP */
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, extractYJ, y, NULL));
    GRB_TRY(GrB_Matrix_clear(delta));
    GRB_TRY(GxB_Matrix_build_Vector(delta, delta_vec,
				    Jvec, delta_vec, GxB_IGNORE_DUP, extract_desc));

    //make delta anti-symmetric: delta_mat = (delta - delta')
    GRB_TRY(GxB_eWiseUnion(delta_mat, NULL, NULL, GrB_MINUS_FP64,   /* FLOP */
			   delta, zero, delta, zero, GrB_DESC_T1));

    //update R
    // R<delta_mat> = UpdateFlows (R, delta_mat) using eWiseMult
    GRB_TRY(GrB_eWiseMult(R, delta_mat, NULL, UpdateFlows,  /* FLOP */
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

  //----------------------------------------------------------------------------
  // optionally construct the output flow matrix, if requested
  //----------------------------------------------------------------------------

  if (flow_mtx != NULL)
  {
    GRB_TRY(GrB_apply(*flow_mtx, NULL, NULL, extractMatrixFlows, R, NULL));
    GRB_TRY(GrB_select(*flow_mtx, NULL, NULL, GrB_VALUEGE_FP64, *flow_mtx, 0, NULL));
  }

  //----------------------------------------------------------------------------
  // for test coverage only
  //----------------------------------------------------------------------------

  #ifdef COVERAGE
  // The MxeAdd operator is not tested via the call to GrB_mxv with the
  // MxeSemiring above, so test it via the MxeAddMonoid.
  GrB_free(&y);
  GRB_TRY(GrB_Vector_new(&y, ResultTuple, 3));
  if (n > NBIG)
  {
    MF_resultTuple64 a = {.d = 1, .j = 2, .residual = 3};
    MF_resultTuple64 b = {.d = 4, .j = 5, .residual = 6};
    GRB_TRY (GrB_Vector_setElement_UDT (y, (void *) &a, 0)) ;
    GRB_TRY (GrB_Vector_setElement_UDT (y, (void *) &b, 0)) ;
    MF_resultTuple64 c = {.d = 0, .j = 0, .residual = 0};
    GRB_TRY (GrB_Vector_reduce_UDT ((void *) &c, NULL, MxeAddMonoid, y, NULL)) ;
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
    GRB_TRY (GrB_Vector_reduce_UDT ((void *) &c, NULL, MxeAddMonoid, y, NULL)) ;
//  printf ("c: resid %g j %g d %g\n", c.residual, (double) c.j, (double) c.d) ;
    LG_ASSERT ((c.residual == 6 && c.j == 5 && c.d == 4), GrB_PANIC) ;
  }
  #endif

  //----------------------------------------------------------------------------
  // free workspace and return result
  //----------------------------------------------------------------------------

  LG_FREE_ALL;
  return GrB_SUCCESS;
#else
  return GrB_NOT_IMPLEMENTED ;
#endif
}
