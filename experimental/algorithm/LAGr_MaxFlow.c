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

// LAGr_MaxFlow is a GraphBLAS implementation of the push-relabel algorithm
// of Baumstark et al. [1]
//
// [1] N. Baumstark, G. E. Blelloch, and J. Shun, "Efficient Implementation of
// a Synchronous Parallel Push-Relabel Algorithm." In: Bansal, N., Finocchi, I.
// (eds) Algorithms - ESA 2015. Lecture Notes in Computer Science(), vol 9294.
// Springer, Berlin, Heidelberg.  https://doi.org/10.1007/978-3-662-48350-3 10.

#include <LAGraphX.h>
#include "LG_internal.h"
#include <LAGraph.h>

#if LG_SUITESPARSE_GRAPHBLAS_V10

//------------------------------------------------------------------------------
// LG_augment_maxflow: sum current excess flow into the output flow
//------------------------------------------------------------------------------

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

    // e<!struct([src,sink])> = select e where (e > 0)
    GRB_TRY(GrB_select(e, src_and_sink, NULL, GrB_VALUEGT_FP64, e, 0, GrB_DESC_RSC));

    // n_active = # of entries in e
    GRB_TRY(GrB_Vector_nvals(n_active, e));
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// LG_global_relabel: global relabeling, based on a BFS from the sink node
//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                        \
{                                           \
    GrB_free(&C);                           \
    GrB_free(&T);                           \
    LAGraph_Delete(&G2, msg);               \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL LG_FREE_WORK

static GrB_Info LG_global_relabel
(
    // inputs:
    GrB_Matrix R,               // flow matrix
    GrB_Index sink,             // sink node
    GrB_Vector src_and_sink,    // mask vector, with just [src sink]
    GrB_UnaryOp GetResidual,    // unary op to compute resid=capacity-flow
    // input/output:
    GrB_Vector d,       // d(i) = height/label of node i
    // outputs:
    GrB_Vector *lvl,    // lvl(i) = distance of node i from sink, if reachable
    char *msg
)
{
    GrB_Matrix C = NULL, T = NULL ;
    LAGraph_Graph G2 = NULL ;
    GrB_Index n ;
    GRB_TRY(GrB_Matrix_nrows(&n, R)) ;
    GRB_TRY(GrB_Matrix_new(&T, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&C, GrB_FP64, n, n));
    // C = GetResidual (R), computing the residual of each edge
    GRB_TRY(GrB_apply(C, NULL, NULL, GetResidual, R, NULL)) ;
    // prune zeros and negative entries from C
    GRB_TRY(GrB_select(C, NULL, NULL, GrB_VALUEGT_FP64, C, 0, NULL)) ;
    // T = C'
    GRB_TRY(GrB_transpose(T, NULL, NULL, C, NULL));
    // construct G2 and its cached transpose and outdegree
    LG_TRY(LAGraph_New(&G2, &T, LAGraph_ADJACENCY_DIRECTED, msg));
    G2->AT = C ;
    C = NULL ;
    LG_TRY(LAGraph_Cached_OutDegree(G2, msg));
    // compute lvl using bfs on G2, starting at sink node
    LG_TRY(LAGr_BreadthFirstSearch(lvl, NULL, G2, sink, msg));
    // d<!struct([src,sink])> = lvl
    GRB_TRY(GrB_assign(d, src_and_sink, NULL, *lvl, GrB_ALL, n, GrB_DESC_SC));
    // d<!struct(lvl)> = n
    GRB_TRY(GrB_assign(d, *lvl, NULL, n, GrB_ALL, n, GrB_DESC_SC));
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------

#undef  LG_FREE_WORK
#define LG_FREE_WORK                        \
{                                           \
    GrB_free(&FlowEdge);                    \
    GrB_free(&CompareTuple);                \
    GrB_free(&ResultTuple);                 \
    GrB_free(&e);                           \
    GrB_free(&d);                           \
    GrB_free(&theta);                       \
    GrB_free(&R);                           \
    GrB_free(&Delta);                       \
    GrB_free(&delta_vec);                   \
    GrB_free(&Map);                         \
    GrB_free(&y);                           \
    GrB_free(&yd);                          \
    GrB_free(&src_and_sink);                \
    GrB_free(&Jvec);                        \
    GrB_free(&Prune);                       \
    GrB_free(&UpdateFlow);                  \
    GrB_free(&Relabel);                     \
    GrB_free(&ResidualFlow);                \
    GrB_free(&MxeIndexMult);                \
    GrB_free(&MxeMult);                     \
    GrB_free(&MxeAdd);                      \
    GrB_free(&MxeAddMonoid);                \
    GrB_free(&MxeSemiring);                 \
    GrB_free(&ExtractJ);                    \
    GrB_free(&CreateCompareVec);            \
    GrB_free(&RxdSemiring);                 \
    GrB_free(&RxdAdd);                      \
    GrB_free(&RxdAddMonoid);                \
    GrB_free(&RxdIndexMult);                \
    GrB_free(&RxdMult);                     \
    GrB_free(&InitForw);                    \
    GrB_free(&InitBack);                    \
    GrB_free(&ResidualForward);             \
    GrB_free(&ResidualBackward);            \
    GrB_free(&zero);                        \
    GrB_free(&empty);                       \
    GrB_free(&t);                           \
    GrB_free(&invariant);                   \
    GrB_free(&CheckInvariant);              \
    GrB_free(&check);                       \
    GrB_free(&ExtractYJ);                   \
    GrB_free(&desc);                        \
    GrB_free(&MakeFlow);                    \
    GrB_free(&GetResidual);                 \
    GrB_free(&lvl) ;                        \
    GrB_free(&ExtractMatrixFlow);           \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL                         \
{                                           \
    LG_FREE_WORK ;                          \
    GrB_free(flow_mtx) ;                    \
}

// helper macro for creating types and operators
#define JIT_STR(f, var) char* var = #f; f

//casting for unary ops
#define F_UNARY(f) ((void (*)(void *, const void *))f)

// casting for index binary ops
#define F_INDEX_BINARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void *, GrB_Index, GrB_Index, const void *)) f)

// casting for binary op
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

//------------------------------------------------------------------------------
// custom types
//------------------------------------------------------------------------------

// type of the R matrix: MF_flowEdge (FlowEdge)
JIT_STR(typedef struct{
  double capacity;      /* original edge weight A(i,j), always positive */
  double flow;          /* current flow along this edge (i,j); can be negative */
  } MF_flowEdge;, FLOWEDGE_STR)

// type of the y vector for y = R*d: MF_resultTuple64/32 (ResultTuple)
JIT_STR(typedef struct{
  double residual;      /* residual = capacity - flow for the edge (i,j) */
  int64_t d;            /* d(j) of the target node j */
  int64_t j;            /* node id of the target node j */
  } MF_resultTuple64;, RESULTTUPLE_STR64)
JIT_STR(typedef struct{
  double residual;      /* residual = capacity - flow for the edge (i,j) */
  int32_t d;            /* d(j) of the target node j */
  int32_t j;            /* node id of the target node j */
  } MF_resultTuple32;, RESULTTUPLE_STR32)

// type of the Map matrix and yd vector: MF_compareTuple64/32 (CompareTuple)
JIT_STR(typedef struct{
  double residual;      /* residual = capacity - flow for the edge (i,j) */
  int64_t di;           /* d(i) for node i */
  int64_t dj;           /* d(j) for node j */
  int64_t j;            /* node id for node j */
  } MF_compareTuple64;, COMPARETUPLE_STR64)
JIT_STR(typedef struct{
  double residual;      /* residual = capacity - flow for the edge (i,j) */
  int32_t di;           /* d(i) for node i */
  int32_t dj;           /* d(j) for node j */
  int32_t j;            /* node id for node j */
  int32_t unused;       /* to pad the struct to 24 bytes */
  } MF_compareTuple32;, COMPARETUPLE_STR32) // 24 bytes: padded

//------------------------------------------------------------------------------
// unary ops to create R from input adjacency matrix G->A and G->AT
//------------------------------------------------------------------------------

// unary op for R = ResidualForward (A)
JIT_STR(void MF_ResidualForward(MF_flowEdge *z, const double *y) {
  z->capacity = (*y);
  z->flow = 0;
  }, CRF_STR)

// unary op for R<!struct(A)> = ResidualBackward (AT)
JIT_STR(void MF_ResidualBackward(MF_flowEdge *z, const double *y) {
  z->capacity = 0;
  z->flow = 0;
  }, CRB_STR)

//------------------------------------------------------------------------------
// R*d semiring
//------------------------------------------------------------------------------

// multiplicative operator, z = R(i,j) * d(j), 64-bit case
JIT_STR(void MF_RxdMult64(MF_resultTuple64 *z,
    const MF_flowEdge *x, GrB_Index i, GrB_Index j,
    const int64_t *y, GrB_Index iy, GrB_Index jy,
    const int64_t* theta) {
  double r = x->capacity - x->flow;
  if(r > 0){
    z->residual = r;
    z->d = (*y);
    z->j = j;
  }
  else{
    z->residual = 0;
    z->d = INT64_MAX;
    z->j = -1;
  }
}, RXDMULT_STR64)

// multiplicative operator, z = R(i,j) * d(j), 32-bit case
JIT_STR(void MF_RxdMult32(MF_resultTuple32 *z,
    const MF_flowEdge *x, GrB_Index i, GrB_Index j,
    const int32_t *y, GrB_Index iy, GrB_Index jy,
    const int32_t* theta) {
  double r = x->capacity - x->flow;
  if(r > 0){
    z->residual = r;
    z->d = (*y);
    z->j = j;
  }
  else{
    z->residual = 0;
    z->d = INT32_MAX;
    z->j = -1;
  }
}, RXDMULT_STR32)

// additive monoid: z = the best tuple, x or y, 64-bit case
JIT_STR(void MF_RxdAdd64(MF_resultTuple64 * z,
    const MF_resultTuple64 * x,
    const MF_resultTuple64 * y) {
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

// additive monoid: z = the best tuple, x or y, 32-bit case
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

//------------------------------------------------------------------------------
// unary ops for delta_vec = ResidualFlow (y)
//------------------------------------------------------------------------------

JIT_STR(void MF_ResidualFlow64(double *z, const MF_resultTuple64 *x)
    { (*z) = x->residual; }, RESIDUALFLOW_STR64)

JIT_STR(void MF_ResidualFlow32(double *z, const MF_resultTuple32 *x)
    { (*z) = x->residual; }, RESIDUALFLOW_STR32)

//------------------------------------------------------------------------------
// binary op for R<Delta> = UpdateFlow (R, Delta) using eWiseMult
//------------------------------------------------------------------------------

JIT_STR(void MF_UpdateFlow(MF_flowEdge *z,
    const MF_flowEdge *x, const double *y) {
  z->capacity = x->capacity;
  z->flow = x->flow + (*y);
  }, UPDATEFLOW_STR)

//------------------------------------------------------------------------------
// binary op for d<struct(y)> = Relabel (d, y) using eWiseMult
//------------------------------------------------------------------------------

JIT_STR(void MF_Relabel64(int64_t *z,
    const int64_t *x, const MF_resultTuple64 *y) {
  if((*x) < y->d+1){
    (*z) = y->d + 1;
  }
  else {
    (*z) = (*x);
  }
  }, RELABEL_STR64)

JIT_STR(void MF_Relabel32(int32_t *z,
    const int32_t *x, const MF_resultTuple32 *y) {
  if((*x) < y->d+1){
    (*z) = y->d + 1;
  }
  else {
    (*z) = (*x);
  }
  }, RELABEL_STR32)

//------------------------------------------------------------------------------
// unary op for Jvec = ExtractJ (yd), where Jvec(i) = yd(i)->j
//------------------------------------------------------------------------------

JIT_STR(void MF_ExtractJ64(int64_t *z, const MF_compareTuple64 *x) { (*z) = x->j; }, EXTRACTJ_STR64)

JIT_STR(void MF_ExtractJ32(int32_t *z, const MF_compareTuple32 *x) { (*z) = x->j; }, EXTRACTJ_STR32)

//------------------------------------------------------------------------------
// unary op for Jvec = ExtractYJ (y), where Jvec(i) = y(i)->j
//------------------------------------------------------------------------------

JIT_STR(void MF_ExtractYJ64(int64_t *z, const MF_resultTuple64 *x) { (*z) = x->j; }, EXTRACTYJ_STR64)

JIT_STR(void MF_ExtractYJ32(int32_t *z, const MF_resultTuple32 *x) { (*z) = x->j; }, EXTRACTYJ_STR32)

//------------------------------------------------------------------------------
// binary op for R(src,:) = InitForw (R (src,:), t')
//------------------------------------------------------------------------------

JIT_STR(void MF_InitForw(MF_flowEdge * z,
    const MF_flowEdge * x, const MF_flowEdge * y){
  z->capacity = x->capacity;
  z->flow = y->flow + x->flow;
  }, INITFORW_STR)

//------------------------------------------------------------------------------
// binary op for R(:,src) = InitBack (R (:,src), t)
//------------------------------------------------------------------------------

JIT_STR(void MF_InitBack(MF_flowEdge * z,
    const MF_flowEdge * x, const MF_flowEdge * y){
  z->capacity = x->capacity;
  z->flow = x->flow - y->flow;
  }, INITBACK_STR)

//------------------------------------------------------------------------------
// y = Map*e semiring
//------------------------------------------------------------------------------

// multiplicative operator, z = Map(i,j)*e(j), 64-bit case
JIT_STR(void MF_MxeMult64(MF_resultTuple64 * z,
    const MF_compareTuple64 * x, GrB_Index i, GrB_Index j,
    const double * y, GrB_Index iy, GrB_Index jy,
    const int64_t* theta){
  bool j_active = ((*y) > 0) ;
  if ((x->di <  x->dj-1) /* case a */
  ||  (x->di == x->dj-1 && !j_active) /* case b */
  ||  (x->di == x->dj   && (!j_active || (j_active && (i < j)))) /* case c */
  ||  (x->di == x->dj+1))   /* case d */
  {
      z->residual = x->residual;
      z->d = x->dj;
      z->j = x->j;
  }
  else
  {
      z->residual = 0;
      z->d = INT64_MAX;
      z->j = -1;
  }
}, MXEMULT_STR64)

// multiplicative operator, z = Map(i,j)*e(j), 32-bit case
JIT_STR(void MF_MxeMult32(MF_resultTuple32 * z,
    const MF_compareTuple32 * x, GrB_Index i, GrB_Index j,
    const double * y, GrB_Index iy, GrB_Index jy,
    const int32_t* theta){
  bool j_active = ((*y) > 0) ;
  if ((x->di <  x->dj-1) /* case a */
  ||  (x->di == x->dj-1 && !j_active) /* case b */
  ||  (x->di == x->dj   && (!j_active || (j_active && (i < j)))) /* case c */
  ||  (x->di == x->dj+1))   /* case d */
  {
      z->residual = x->residual;
      z->d = x->dj;
      z->j = x->j;
  }
  else
  {
      z->residual = 0;
      z->d = INT32_MAX;
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

//------------------------------------------------------------------------------
// binary op for yd = CreateCompareVec (y,d) using eWiseMult
//------------------------------------------------------------------------------

JIT_STR(void MF_CreateCompareVec64(MF_compareTuple64 *comp,
    const MF_resultTuple64 *res, const int64_t *height) {
  comp->di = (*height);
  comp->residual = res->residual;
  comp->dj = res->d;
  comp->j = res->j;
  }, CREATECOMPAREVEC_STR64)

JIT_STR(void MF_CreateCompareVec32(MF_compareTuple32 *comp,
    const MF_resultTuple32 *res, const int32_t *height) {
  comp->di = (*height);
  comp->residual = res->residual;
  comp->dj = res->d;
  comp->j = res->j;
  comp->unused = 0 ;
  }, CREATECOMPAREVEC_STR32)

//------------------------------------------------------------------------------
// index unary op to remove empty tuples from y (for which y->j is -1)
//------------------------------------------------------------------------------

JIT_STR(void MF_Prune64(bool * z, const MF_resultTuple64 * x,
  GrB_Index ix, GrB_Index jx, const int64_t * theta){
  *z = (x->j != -1) ;
  }, PRUNE_STR64)

JIT_STR(void MF_Prune32(bool * z, const MF_resultTuple32 * x,
  GrB_Index ix, GrB_Index jx, const int32_t * theta){
  *z = (x->j != -1) ;
  }, PRUNE_STR32)

//------------------------------------------------------------------------------
// unary op for t = MakeFlow (e), where t(i) = (0, e(i))
//------------------------------------------------------------------------------

JIT_STR(void MF_MakeFlow(MF_flowEdge * flow_edge, const double * flow){
  flow_edge->capacity = 0;
  flow_edge->flow = (*flow);
  }, MAKEFLOW_STR)

//------------------------------------------------------------------------------
// binary op CheckInvariant to check invariants (debugging only)
//------------------------------------------------------------------------------

#ifdef DBG
JIT_STR(void MF_CheckInvariant64(bool *z, const int64_t *height,
    const MF_resultTuple64 *result) {
  (*z) = ((*height) == result->d+1);
  }, CHECKINVARIANT_STR64)

JIT_STR(void MF_CheckInvariant32(bool *z, const int32_t *height,
    const MF_resultTuple32 *result) {
  (*z) = ((*height) == result->d+1);
  }, CHECKINVARIANT_STR32)
#endif

//------------------------------------------------------------------------------
// binary op for C = GetResidual (R), computing the residual of each edge
//------------------------------------------------------------------------------

JIT_STR(void MF_GetResidual(double * res, const MF_flowEdge * flow_edge){
    (*res) = flow_edge->capacity - flow_edge->flow;
}, GETRESIDUAL_STR)

//------------------------------------------------------------------------------
// unary op for flow_mtx = ExtractMatrixFlow (R)
//------------------------------------------------------------------------------

JIT_STR(void MF_ExtractMatrixFlow(double* flow, const MF_flowEdge* edge){*flow = edge->flow;}, EMFLOW_STR)

#endif

//------------------------------------------------------------------------------
// LAGraph_MaxFlow
//------------------------------------------------------------------------------

int LAGr_MaxFlow
(
    // output:
    double *f,              // max flow from src node to sink node
    GrB_Matrix *flow_mtx,   // optional output flow matrix
    // input:
    LAGraph_Graph G,        // graph to compute maxflow on
    GrB_Index src,          // source node
    GrB_Index sink,         // sink node
    char *msg
)
{

#if LG_SUITESPARSE_GRAPHBLAS_V10

  //----------------------------------------------------------------------------
  // declare variables
  //----------------------------------------------------------------------------

  // types
  GrB_Type FlowEdge = NULL ;
  GrB_Type ResultTuple = NULL ;
  GrB_Type CompareTuple = NULL ;

  GrB_Vector lvl = NULL ;
  GrB_UnaryOp GetResidual = NULL ;

  // to create R
  GrB_UnaryOp ResidualForward = NULL, ResidualBackward = NULL ;
  GrB_Matrix R = NULL ;

  // to initialize R with initial saturated flows
  GrB_Vector e = NULL, t = NULL ;
  GrB_UnaryOp MakeFlow = NULL ;
  GrB_BinaryOp InitForw = NULL, InitBack = NULL ;

  // height/label vector
  GrB_Vector d = NULL ;

  // src and sink mask vector and n_active
  GrB_Vector src_and_sink = NULL ;
  GrB_Index n_active = INT64_MAX ;

  // semiring and vectors for y<struct(e)> = R x d
  GrB_Vector y = NULL ;
  GrB_IndexUnaryOp Prune = NULL ;
  GxB_IndexBinaryOp RxdIndexMult = NULL ;
  GrB_BinaryOp RxdAdd = NULL, RxdMult = NULL ;
  GrB_Monoid RxdAddMonoid = NULL ;
  GrB_Semiring RxdSemiring = NULL ;
  GrB_Scalar theta = NULL ;

  // binary op and yd
  GrB_Vector yd = NULL ;
  GrB_BinaryOp CreateCompareVec = NULL ;

  // utility vectors, Matrix, and ops for mapping
  GrB_Matrix Map = NULL ;
  GrB_Vector Jvec = NULL ;
  GrB_UnaryOp ExtractJ = NULL, ExtractYJ = NULL ;

  // Map*e semiring
  GrB_Semiring MxeSemiring = NULL ;
  GrB_Monoid MxeAddMonoid = NULL ;
  GrB_BinaryOp MxeAdd = NULL, MxeMult = NULL ;
  GxB_IndexBinaryOp MxeIndexMult = NULL ;

  // to extract the residual flow
  GrB_UnaryOp ResidualFlow = NULL ;
  GrB_UnaryOp ExtractMatrixFlow = NULL ;

  // Delta structures
  GrB_Vector delta_vec = NULL ;
  GrB_Matrix Delta = NULL ;

  // update height
  GrB_BinaryOp Relabel = NULL ;

  // update R structure
  GrB_BinaryOp UpdateFlow = NULL ;

  // scalars
  GrB_Scalar zero = NULL ;
  GrB_Scalar empty = NULL ;

  // invariant (for debugging only)
  GrB_Vector invariant = NULL ;
  GrB_BinaryOp CheckInvariant = NULL ;
  GrB_Scalar check = NULL ;
  bool check_raw;

  // descriptor for matrix building
  GrB_Descriptor desc = NULL ;

  //----------------------------------------------------------------------------
  // check inputs
  //----------------------------------------------------------------------------

  if (flow_mtx != NULL)
  {
    (*flow_mtx) = NULL ;
  }
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

  // create types for computation
  GRB_TRY(GxB_Type_new(&FlowEdge, sizeof(MF_flowEdge), "MF_flowEdge", FLOWEDGE_STR));

  GRB_TRY(GxB_UnaryOp_new(&GetResidual, F_UNARY(MF_GetResidual), GrB_FP64, FlowEdge,
        "MF_GetResidual", GETRESIDUAL_STR));

  #ifdef DBG
  GRB_TRY(GrB_Scalar_new(&check, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement(check, false));
  GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
  #endif

  // ops create R from A
  GRB_TRY(GxB_UnaryOp_new(&ResidualForward,
        F_UNARY(MF_ResidualForward), FlowEdge , GrB_FP64,
        "MF_ResidualForward", CRF_STR));
  GRB_TRY(GxB_UnaryOp_new(&ResidualBackward,
        F_UNARY(MF_ResidualBackward), FlowEdge , GrB_FP64,
        "MF_ResidualBackward", CRB_STR));

  // ops to initialize R with initial saturated flows from the source node
  GRB_TRY(GxB_BinaryOp_new(&InitForw,
        F_BINARY(MF_InitForw), FlowEdge, FlowEdge, FlowEdge,
        "MF_InitForw", INITFORW_STR));
  GRB_TRY(GxB_BinaryOp_new(&InitBack,
        F_BINARY(MF_InitBack), FlowEdge, FlowEdge, FlowEdge,
        "MF_InitBack", INITBACK_STR));
  GRB_TRY(GxB_UnaryOp_new(&MakeFlow, F_UNARY(MF_MakeFlow), FlowEdge, GrB_FP64,
        "MF_MakeFlow", MAKEFLOW_STR));

  // construct [src,sink] mask
  GRB_TRY(GrB_Vector_new(&src_and_sink, GrB_BOOL, n));
  GRB_TRY (GrB_Vector_setElement (src_and_sink, true, sink)) ;
  GRB_TRY (GrB_Vector_setElement (src_and_sink, true, src)) ;

  // create delta vector and Delta matrix
  GRB_TRY(GrB_Matrix_new(&Delta, GrB_FP64, n, n));
  GRB_TRY(GrB_Vector_new(&delta_vec, GrB_FP64, n));

  // operator to update R structure
  GRB_TRY(GxB_BinaryOp_new(&UpdateFlow,
        F_BINARY(MF_UpdateFlow), FlowEdge, FlowEdge, GrB_FP64,
        "MF_UpdateFlow", UPDATEFLOW_STR));

  // create scalars
  GRB_TRY(GrB_Scalar_new(&zero, GrB_FP64));
  GRB_TRY(GrB_Scalar_setElement(zero, 0));
  GRB_TRY(GrB_Scalar_new (&empty, GrB_FP64)) ;

  // create op for optional output flow_mtx
  if (flow_mtx != NULL)
  {
    GRB_TRY(GxB_UnaryOp_new(&ExtractMatrixFlow,
        F_UNARY(MF_ExtractMatrixFlow), GrB_FP64, FlowEdge,
        "MF_ExtractMatrixFlow", EMFLOW_STR));
  }

  //----------------------------------------------------------------------------
  // determine the integer type to use for the problem
  //----------------------------------------------------------------------------

  GrB_Type Integer_Type = NULL ;

  #ifdef COVERAGE
  // Just for test coverage, use 64-bit ints for n > 100.  Do not use this
  // rule in production!
  #define NBIG 100
  #else
  // For production use: 64-bit integers if n > 2^31
  #define NBIG INT32_MAX
  #endif
  if (n > NBIG){

    //--------------------------------------------------------------------------
    // use 64-bit integers
    //--------------------------------------------------------------------------

    Integer_Type = GrB_INT64 ;

    // create types for computation
    GRB_TRY(GxB_Type_new(&ResultTuple, sizeof(MF_resultTuple64),
        "MF_resultTuple64", RESULTTUPLE_STR64));
    GRB_TRY(GxB_Type_new(&CompareTuple, sizeof(MF_compareTuple64),
        "MF_compareTuple64", COMPARETUPLE_STR64));

    // invariant check
    #ifdef DBG
    GRB_TRY(GxB_BinaryOp_new(&CheckInvariant,
        F_BINARY(MF_CheckInvariant64), GrB_BOOL, GrB_INT64, ResultTuple,
        "MF_CheckInvariant64", CHECKINVARIANT_STR64));
    #endif

    GRB_TRY(GxB_UnaryOp_new(&ResidualFlow,
        F_UNARY(MF_ResidualFlow64), GrB_FP64, ResultTuple,
        "MF_ResidualFlow64", RESIDUALFLOW_STR64));

    // create ops for R*d semiring
    GRB_TRY(GrB_Scalar_new(&theta, GrB_INT64));
    GRB_TRY(GrB_Scalar_setElement_INT64(theta, 0));
    GRB_TRY(GxB_IndexBinaryOp_new(&RxdIndexMult,
        F_INDEX_BINARY(MF_RxdMult64), ResultTuple, FlowEdge, GrB_INT64, GrB_INT64,
        "MF_RxdMult64", RXDMULT_STR64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&RxdMult, RxdIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&RxdAdd,
        F_BINARY(MF_RxdAdd64), ResultTuple, ResultTuple, ResultTuple,
        "MF_RxdAdd64", RXDADD_STR64));
    MF_resultTuple64 id = {.d = INT64_MAX, .j = -1, .residual = 0};
    GRB_TRY(GrB_Monoid_new_UDT(&RxdAddMonoid, RxdAdd, &id));

    // create binary op for yd
    GRB_TRY(GxB_BinaryOp_new(&CreateCompareVec,
        F_BINARY(MF_CreateCompareVec64), CompareTuple, ResultTuple, GrB_INT64,
        "MF_CreateCompareVec64", CREATECOMPAREVEC_STR64));

    // create op to prune empty tuples
    GRB_TRY(GxB_IndexUnaryOp_new(&Prune,
        (GxB_index_unary_function) MF_Prune64, GrB_BOOL, ResultTuple, GrB_INT64,
        "MF_Prune64", PRUNE_STR64));

    // create ops for mapping
    GRB_TRY(GxB_UnaryOp_new(&ExtractJ,
        F_UNARY(MF_ExtractJ64), GrB_INT64, CompareTuple,
        "MF_ExtractJ64", EXTRACTJ_STR64));
    GRB_TRY(GxB_UnaryOp_new(&ExtractYJ,
        F_UNARY(MF_ExtractYJ64), GrB_INT64, ResultTuple,
        "MF_ExtractYJ64", EXTRACTYJ_STR64));

    // create ops for Map*e semiring
    GRB_TRY(GxB_IndexBinaryOp_new(&MxeIndexMult,
        F_INDEX_BINARY(MF_MxeMult64), ResultTuple, CompareTuple, GrB_FP64, GrB_INT64,
        "MF_MxeMult64", MXEMULT_STR64));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MxeMult, MxeIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&MxeAdd,
        F_BINARY(MF_MxeAdd64), ResultTuple, ResultTuple, ResultTuple,
        "MF_MxeAdd64", MXEADD_STR64));
    GRB_TRY(GrB_Monoid_new_UDT(&MxeAddMonoid, MxeAdd, &id));

    // update height binary op
    GRB_TRY(GxB_BinaryOp_new(&Relabel,
        F_BINARY(MF_Relabel64), GrB_INT64, GrB_INT64, ResultTuple,
        "MF_Relabel64", RELABEL_STR64));

  }else{

    //--------------------------------------------------------------------------
    // use 32-bit integers
    //--------------------------------------------------------------------------

    Integer_Type = GrB_INT32 ;

    // create types for computation
    GRB_TRY(GxB_Type_new(&ResultTuple, sizeof(MF_resultTuple32),
        "MF_resultTuple32", RESULTTUPLE_STR32));
    GRB_TRY(GxB_Type_new(&CompareTuple, sizeof(MF_compareTuple32),
        "MF_compareTuple32", COMPARETUPLE_STR32));

    // invariant check
    #ifdef DBG
    GRB_TRY(GxB_BinaryOp_new(&CheckInvariant,
        F_BINARY(MF_CheckInvariant32), GrB_BOOL, GrB_INT32, ResultTuple,
        "MF_CheckInvariant32", CHECKINVARIANT_STR32));
    #endif

    GRB_TRY(GxB_UnaryOp_new(&ResidualFlow,
        F_UNARY(MF_ResidualFlow32), GrB_FP64, ResultTuple,
        "MF_ResidualFlow32", RESIDUALFLOW_STR32));

    // create ops for R*d semiring
    GRB_TRY(GrB_Scalar_new(&theta, GrB_INT32));
    GRB_TRY(GrB_Scalar_setElement_INT32(theta, 0));
    GRB_TRY(GxB_IndexBinaryOp_new(&RxdIndexMult,
        F_INDEX_BINARY(MF_RxdMult32), ResultTuple, FlowEdge, GrB_INT32, GrB_INT32,
        "MF_RxdMult32", RXDMULT_STR32));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&RxdMult, RxdIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&RxdAdd,
        F_BINARY(MF_RxdAdd32), ResultTuple, ResultTuple, ResultTuple,
        "MF_RxdAdd32", RXDADD_STR32));
    MF_resultTuple32 id = {.d = INT32_MAX, .j = -1, .residual = 0};
    GRB_TRY(GrB_Monoid_new_UDT(&RxdAddMonoid, RxdAdd, &id));

    // create binary op for yd
    GRB_TRY(GxB_BinaryOp_new(&CreateCompareVec,
        F_BINARY(MF_CreateCompareVec32), CompareTuple, ResultTuple, GrB_INT32,
        "MF_CreateCompareVec32", CREATECOMPAREVEC_STR32));

    // create op to prune empty tuples
    GRB_TRY(GxB_IndexUnaryOp_new(&Prune,
        (GxB_index_unary_function) MF_Prune32, GrB_BOOL, ResultTuple, GrB_INT32,
        "MF_Prune32", PRUNE_STR32));

    // create ops for mapping
    GRB_TRY(GxB_UnaryOp_new(&ExtractJ,
        F_UNARY(MF_ExtractJ32), GrB_INT32, CompareTuple,
        "MF_ExtractJ32", EXTRACTJ_STR32));
    GRB_TRY(GxB_UnaryOp_new(&ExtractYJ,
        F_UNARY(MF_ExtractYJ32), GrB_INT32, ResultTuple,
        "MF_ExtractYJ32", EXTRACTYJ_STR32));

    // create ops for Map*e semiring
    GRB_TRY(GxB_IndexBinaryOp_new(&MxeIndexMult,
        F_INDEX_BINARY(MF_MxeMult32), ResultTuple, CompareTuple, GrB_FP64, GrB_INT32,
        "MF_MxeMult32", MXEMULT_STR32));
    GRB_TRY(GxB_BinaryOp_new_IndexOp(&MxeMult, MxeIndexMult, theta));
    GRB_TRY(GxB_BinaryOp_new(&MxeAdd,
        F_BINARY(MF_MxeAdd32), ResultTuple, ResultTuple, ResultTuple,
        "MF_MxeAdd32", MXEADD_STR32));
    GRB_TRY(GrB_Monoid_new_UDT(&MxeAddMonoid, MxeAdd, &id));

    // update height binary op
    GRB_TRY(GxB_BinaryOp_new(&Relabel,
        F_BINARY(MF_Relabel32), GrB_INT32, GrB_INT32, ResultTuple,
        "MF_Relabel32", RELABEL_STR32));
  }

  //----------------------------------------------------------------------------
  // create remaining vectors, matrices, descriptor, and semirings
  //----------------------------------------------------------------------------

  GRB_TRY(GrB_Matrix_new(&Map, CompareTuple, n,n));
  GRB_TRY(GrB_Vector_new(&Jvec, Integer_Type, n));
  GRB_TRY(GrB_Vector_new(&yd, CompareTuple, n));
  GRB_TRY(GrB_Vector_new(&y, ResultTuple, n));

  GRB_TRY(GrB_Semiring_new(&RxdSemiring, RxdAddMonoid, RxdMult));
  GRB_TRY(GrB_Semiring_new(&MxeSemiring, MxeAddMonoid, MxeMult));

  // create descriptor for building the Map and Delta matrices
  GRB_TRY(GrB_Descriptor_new(&desc));
  GRB_TRY(GrB_set(desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST));

  // create and init d vector
  GRB_TRY(GrB_Vector_new(&d, Integer_Type, n));
  GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(d, NULL, NULL, n, &src, 1, NULL));

  // create R, with no flow
  GRB_TRY(GrB_Matrix_new(&R, FlowEdge, n, n));
  // R = ResidualForward (A)
  GRB_TRY(GrB_apply(R, NULL, NULL, ResidualForward, A, NULL));
  // R<!struct(A)> = ResidualBackward (AT)
  GRB_TRY(GrB_apply(R, A, NULL, ResidualBackward, AT, GrB_DESC_SC));

  // initial global relabeling
  LG_TRY (LG_global_relabel (R, sink, src_and_sink, GetResidual, d, &lvl, msg)) ;

  // create excess vector e and initial flows from the src to its neighbors
  // e<struct(lvl)> = A (src,:)
  GRB_TRY(GrB_Vector_new(&e, GrB_FP64, n));
  GRB_TRY(GrB_extract(e, lvl, NULL, A, GrB_ALL, n, src, GrB_DESC_ST0));
  GrB_free(&lvl);
  // t = MakeFlow (e), where t(i) = (0, e(i))
  GRB_TRY(GrB_Vector_new(&t, FlowEdge, n));
  GRB_TRY(GrB_apply(t, NULL, NULL, MakeFlow, e, NULL));
  // R(src,:) = InitForw (R (src,:), t')
  GRB_TRY(GrB_assign(R, NULL, InitForw, t, src, GrB_ALL, n, NULL));
  // R(:,src) = InitBack (R (:,src), t)
  GRB_TRY(GrB_assign(R, NULL, InitBack, t, GrB_ALL, n, src, NULL));
  GrB_free(&t) ;

  // augment the maxflow with the initial flows from the src to its neighbors
  LG_TRY (LG_augment_maxflow (f, e, sink, src_and_sink, &n_active, msg)) ;

  //----------------------------------------------------------------------------
  // compute the max flow
  //----------------------------------------------------------------------------

  for (int64_t iter = 0 ; n_active > 0 ; iter++)
  {

    //--------------------------------------------------------------------------
    // Part 1: global relabeling
    //--------------------------------------------------------------------------

    if ((iter > 0) && (flow_mtx != NULL) && (iter % 12 == 0))
    {
      LG_TRY (LG_global_relabel (R, sink, src_and_sink, GetResidual, d, &lvl, msg)) ;
      // delete nodes in e that cannot be reached from the sink
      // e<!struct(lvl)> = empty scalar
      GrB_assign (e, lvl, NULL, empty, GrB_ALL, n, GrB_DESC_SC) ;
      GrB_free(&lvl);
      GRB_TRY(GrB_Vector_nvals(&n_active, e));
      if(n_active == 0) break;
    }

    //--------------------------------------------------------------------------
    // Part 2: deciding where to push
    //--------------------------------------------------------------------------

    // y<struct(e),replace> = R*d using the RxdSemiring
    GRB_TRY(GrB_mxv(y, e, NULL, RxdSemiring, R, d, GrB_DESC_RS));

    // remove empty tuples (0,inf,-1) from y
    GRB_TRY(GrB_select(y, NULL, NULL, Prune, y, 0, NULL));

    //--------------------------------------------------------------------------
    // Part 3: verifying the pushes
    //--------------------------------------------------------------------------

    // create Map matrix from pattern and values of yd
    // yd = CreateCompareVec (y,d) using eWiseMult
    GRB_TRY(GrB_eWiseMult(yd, NULL, NULL, CreateCompareVec, y,  d, NULL));
    // Jvec = ExtractJ (yd), where Jvec(i) = yd(i)->j
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, ExtractJ, yd, NULL));
    GRB_TRY(GrB_Matrix_clear(Map));
    GRB_TRY(GrB_Matrix_build(Map, yd, Jvec, yd, GxB_IGNORE_DUP, desc));

    // make e dense for Map computation
    // TODO: consider keeping e in bitmap/full format only,
    // or always full with e(i)=0 denoting a non-active node.
    GRB_TRY(GrB_assign(e, e, NULL, 0, GrB_ALL, n, GrB_DESC_SC));

    // y = Map*e using the MxeSemiring
    GRB_TRY(GrB_mxv(y, NULL, NULL, MxeSemiring, Map, e, NULL));

    // remove empty tuples (0,inf,-1) from y
    GRB_TRY(GrB_select(y, NULL, NULL, Prune, y, -1, NULL));

    // relabel, updating the height/label vector d
    // d<struct(y)> = Relabel (d, y) using eWiseMult
    GRB_TRY(GrB_eWiseMult(d, y, NULL, Relabel, d, y, GrB_DESC_S));

    #ifdef DBG
        // assert invariant for all labels
        GRB_TRY(GrB_eWiseMult(invariant, y, NULL, CheckInvariant, d, y, GrB_DESC_RS));
        GRB_TRY(GrB_reduce(check, NULL, GrB_LAND_MONOID_BOOL, invariant, NULL));
        GRB_TRY(GrB_Scalar_extractElement(&check_raw, check));
        ASSERT(check_raw == true);
    #endif

    //--------------------------------------------------------------------------
    // Part 4: executing the pushes
    //--------------------------------------------------------------------------

    // extract residual flows from y
    // delta_vec = ResidualFlow (y), obtaining just the residual flows
    GRB_TRY(GrB_apply(delta_vec, NULL, NULL, ResidualFlow, y, NULL));

    // delta_vec = min (delta_vec, e), where e is dense
    GRB_TRY(GrB_eWiseMult(delta_vec, NULL, NULL, GrB_MIN_FP64, delta_vec, e, NULL));

    // create the Delta matrix from delta_vec and y
    // note that delta_vec has the same structure as y
    // Jvec = ExtractYJ (y), where Jvec(i) = y(i)->j
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, ExtractYJ, y, NULL));
    GRB_TRY(GrB_Matrix_clear(Delta));
    GRB_TRY(GrB_Matrix_build(Delta, delta_vec, Jvec, delta_vec, GxB_IGNORE_DUP, desc));

    // make Delta anti-symmetric
    // Delta = (Delta - Delta')
    GRB_TRY(GxB_eWiseUnion(Delta, NULL, NULL, GrB_MINUS_FP64, Delta, zero, Delta, zero, GrB_DESC_T1));

    // update R
    // R<Delta> = UpdateFlow (R, Delta) using eWiseMult
    GRB_TRY(GrB_eWiseMult(R, Delta, NULL, UpdateFlow, R, Delta, GrB_DESC_S));

    // reduce Delta to delta_vec
    // delta_vec = sum (Delta), summing up each row of Delta
    GRB_TRY(GrB_reduce(delta_vec, NULL, NULL, GrB_PLUS_FP64, Delta, GrB_DESC_T0));

    // add delta_vec to e
    // e<struct(delta_vec)> += delta_vec
    GRB_TRY(GrB_assign(e, delta_vec, GrB_PLUS_FP64, delta_vec, GrB_ALL, n, GrB_DESC_S));

    // augment maxflow for all active nodes
    LG_TRY (LG_augment_maxflow (f, e, sink, src_and_sink, &n_active, msg)) ;
  }

  //----------------------------------------------------------------------------
  // optionally construct the output flow matrix, if requested
  //----------------------------------------------------------------------------

  if (flow_mtx != NULL)
  {
    // flow_mtx = ExtractMatrixFlow (R)
    GRB_TRY(GrB_Matrix_new(flow_mtx, GrB_FP64, n, n));
    GRB_TRY(GrB_apply(*flow_mtx, NULL, NULL, ExtractMatrixFlow, R, NULL));
    // delete any zero or negative flows from the flow_mtx
    GRB_TRY(GrB_select(*flow_mtx, NULL, NULL, GrB_VALUEGT_FP64, *flow_mtx, 0, NULL));
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
    LG_ASSERT ((c.residual == 6 && c.j == 5 && c.d == 4), GrB_PANIC) ;
  }
  #endif

  //----------------------------------------------------------------------------
  // free workspace and return result
  //----------------------------------------------------------------------------

  LG_FREE_WORK ;
  return GrB_SUCCESS;
#else
  return GrB_NOT_IMPLEMENTED ;
#endif
}
