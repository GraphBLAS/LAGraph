
#include <LG_internal.h>
#include <LAGraphX.h>
#include <LAGraph.h>
#include <math.h>

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK                        \
{                                           \
  GrB_Type_free(&GrB_ResidualEdge);         \
}                                           \

#define LG_FREE_ALL                                 \
{                                                   \
    LG_FREE_WORK ;                                  \
}

//casting for unary ops
#define F_UNARY(f) ((void (*)(void *, const void *))f)

// casting for index binary ops
#define F_INDEX_BINARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void *, GrB_Index, GrB_Index, const void *)) f)

// casting for index unary op
/* #define F_INDEX_UNARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void*)) f) */

// casting for binary op
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

// strings for JIT
#define GRB_FLOWEDGE_STR "typedef struct{ float flow; float capacity; } GrB_Flow_Edge;"
#define GRB_MFRESULT_STR "typedef struct{ float residual; int d; GrB_Index j; } MF_result_tuple; "
#define GRB_CR_STR "void CreateResidual_UOp(GrB_Flow_Edge *f, const float *cap) {f->flow = 0; f->capacity = (*cap);}"

#define GZB_MULT_STR "void Rxd_MultBOp(MF_result_tuple *y, const GrB_Flow_Edge *R, GrB_Index ix, GrB_Index jx, const int *d, GrB_Index iy, GrB_Index ik, const int* theta) { float r = R->capacity - R->flow; if(r > 0){ y->residual = r; y->d = *d; y->j = jx; } else{ y->residual = 0; y->d = INT32_MAX; y->j = -1; } }"
#define GRB_ADD_STR "void Rxd_AddMonoid(MF_result_tuple * z, const MF_result_tuple * x, const MF_result_tuple * y) {if(x->d < y->d){ memcpy(z, x, sizeof(MF_result_tuple)); } else if(x->d > y->d){ memcpy(z, y, sizeof(MF_result_tuple));}else{if(x->residual > y->residual){ memcpy(z, x, sizeof(MF_result_tuple)); } else if(x->residual < y->residual){ memcpy(z, y, sizeof(MF_result_tuple)); } else{ if(x->j > y->j){ memcpy(z, x, sizeof(MF_result_tuple)); } else{ memcpy(z, y, sizeof(MF_result_tuple)); } } } }"
#define GRB_INIT_FLOW_STR "void GrB_init_flows(GrB_Flow_Edge * z, const GrB_Flow_Edge * y, const float * x){z->flow = *x;}"
#define GRB_PRUNE_STR "void prune(int *z, const MF_result_tuple *x) {if(x->d < INT32_MAX){(*z) = x->j;}}"
#define GRB_MIN_STR "void MF_eWiseMin(float *z, const MF_result_tuple *y, const float *x) {if(y->d < INT32_MAX){(*z) = fmin(y->residual, (*x));}}"
#define GRB_UPDATE_STR "void MF_updateFlow(GrB_Flow_Edge *z, const GrB_Flow_Edge *y, const float *x) {(*z) = y->flow + (*x);}"
#define MF_UPDATE_HEIGHT_STR "void MF_updateHeight(int *z, const int *x, const MF_result_tuple *y) {(*z) = y->d;}"

//custom types
typedef struct{
  float flow;
  float capacity;
} GrB_Flow_Edge;

typedef struct{
  float residual;
  int d;
  GrB_Index j;
} MF_result_tuple;

//Make this an Unary apply op
//takes capacity from adj and creates flows for the forward edges 
void CreateResidual_UOp(GrB_Flow_Edge *f, const float *cap) {
  f->flow = 0;
  f->capacity = (*cap);
}

void Rxd_MultBOp(MF_result_tuple *y, const GrB_Flow_Edge *R, GrB_Index ix, GrB_Index jx, const int *d, GrB_Index iy, GrB_Index jy, const int* theta) {
  float r = R->capacity - R->flow;
  if(r > 0){
    y->residual = r;
    y->d = *d;
    y->j = jx;
  }
  else{
    y->residual = 0;
    y->d = INT32_MAX;
    y->j = -1;
  }
}

void MF_Extract(float *z, const float *y, const float *x) {
  if(*x > 0){
    *z = *x;
  }
}

void MF_updateFlow(GrB_Flow_Edge *z, const GrB_Flow_Edge *y, const float *x) {
  z->flow = y->flow + (*x);
}

void MF_updateHeight(int *z, const int *x, const MF_result_tuple *y) {
  (*z) = y->d;
}

void Rxd_AddMonoid(MF_result_tuple * z, const MF_result_tuple * x, const MF_result_tuple * y) {
  if(x->d < y->d){
    memcpy(z, x, sizeof(MF_result_tuple));
  }
  else if(x->d > y->d){
    memcpy(z, y, sizeof(MF_result_tuple));
  }
  else{
    if(x->residual > y->residual){
      memcpy(z, x, sizeof(MF_result_tuple));
    }
    else if(x->residual < y->residual){
      memcpy(z, y, sizeof(MF_result_tuple));
    }
    else{
      if(x->j > y->j){
	memcpy(z, x, sizeof(MF_result_tuple));
      }
      else{
	memcpy(z, y, sizeof(MF_result_tuple));
      }
    }
  }
}

void prune(int *z, const MF_result_tuple *x) {
  if(x->d < INT32_MAX){
    (*z) = x->j;
  }
}

void MF_eWiseMin(float *z, const MF_result_tuple *y, const float *x) {
  if(y->d < INT32_MAX){
    (*z) = fmin(y->residual, (*x));
  }
}

void GrB_init_flows(GrB_Flow_Edge * z, const GrB_Flow_Edge * y, const float * x){
  z->flow = *x;
}
// R is resulting residual graph
// f is max flow
int LAGraph_MaxFlow(LAGraph_Graph G, GrB_Index S, GrB_Index T, int * f, char *msg){
  
  //create semirings and types
  GrB_Type GrB_ResidualEdge;
  GRB_TRY(GxB_Type_new(&GrB_ResidualEdge, sizeof(GrB_Flow_Edge), "GrB_Flow_Edge", GRB_FLOWEDGE_STR));

  GrB_Type GrB_MFResult;
  GRB_TRY(GxB_Type_new(&GrB_MFResult, sizeof(MF_result_tuple), "MF_result_tuple", GRB_MFRESULT_STR));

  GrB_UnaryOp GrB_CR_UOp;
  GRB_TRY(GxB_UnaryOp_new(&GrB_CR_UOp, F_UNARY(CreateResidual_UOp), GrB_ResidualEdge, GrB_FP32, "CreateResidual_UOp", GRB_CR_STR));

  GrB_UnaryOp GrB_prune;
  GRB_TRY(GxB_UnaryOp_new(&GrB_prune, F_UNARY(prune), GrB_INT32, GrB_MFResult, "prune", GRB_PRUNE_STR));

  GrB_BinaryOp init_flow;
  GRB_TRY(GxB_BinaryOp_new(&init_flow, F_BINARY(GrB_init_flows), GrB_ResidualEdge, GrB_ResidualEdge, GrB_FP32, "GrB_init_flows", GRB_INIT_FLOW_STR)); //accum

  GrB_BinaryOp GrB_MFeWiseMin;
  GRB_TRY(GxB_BinaryOp_new(&GrB_MFeWiseMin, F_BINARY(MF_eWiseMin), GrB_FP32, GrB_MFResult, GrB_FP32, "MF_eWiseMin", GRB_MIN_STR));

  GzB_IndexBinaryOp GzB_MFMult;
  GrB_BinaryOp GrB_RxdMult;
  GrB_Scalar theta;
  GRB_TRY(GrB_Scalar_new(&theta, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement_INT32(theta, 1));
  GRB_TRY(GzB_IndexBinaryOp_new(&GzB_MFMult, F_INDEX_BINARY(Rxd_MultBOp), GrB_MFResult, GrB_ResidualEdge, GrB_INT32, GrB_INT32,  "Rxd_MultBOp", GZB_MULT_STR));
  GRB_TRY(GzB_BinaryOp_new_IndexOp(&GrB_RxdMult, GzB_MFMult, theta));

  GrB_BinaryOp GrB_MFAdd;
  GrB_Monoid GrB_RxdAdd;
  GRB_TRY(GxB_BinaryOp_new(&GrB_MFAdd, F_BINARY(Rxd_AddMonoid), GrB_MFResult, GrB_MFResult, GrB_MFResult, "Rxd_AddMonoid", GRB_ADD_STR));
  MF_result_tuple id;
  id.d = INT32_MAX;
  id.j = -1;
  id.residual = 0;
  GRB_TRY(GrB_Monoid_new(&GrB_RxdAdd, GrB_MFAdd, (void*)&id));

  GrB_Semiring Rxd_semiring;
  GRB_TRY(GrB_Semiring_new(&Rxd_semiring, GrB_RxdAdd, GrB_RxdMult));

  //make R symmetyric of resiual edge type
  GrB_Matrix A = G->A;
  GrB_Matrix R = NULL;
  GrB_Index n;
  GRB_TRY(GrB_Matrix_nrows(&n, A));
  GRB_TRY(GrB_Matrix_new(&R, GrB_ResidualEdge, n, n));
  GRB_TRY(GrB_apply(R, NULL, NULL, GrB_CR_UOp, A, NULL));
  GRB_TRY(GrB_apply(R, R, NULL, GrB_CR_UOp, A, GrB_DESC_SCT1));
  
  //create d (height) vector and e (excess) vector
  GrB_Vector d = NULL;
  GrB_Vector e = NULL;
  GRB_TRY(GrB_Vector_new(&e, GrB_FP32, n));
  GRB_TRY(GrB_Vector_new(&d, GrB_INT32, n));
  GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));

  //init e and d
  GrB_Vector t;
  GRB_TRY(GrB_Vector_new(&t, GrB_FP32, n));
  GrB_Scalar size;
  GRB_TRY(GrB_Scalar_new(&size, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement_INT32(size, n));
  GRB_TRY(GrB_Vector_setElement(d, size, S));
  GRB_TRY(GrB_Vector_setElement(t, 1, S)); 
  GRB_TRY(GrB_mxv(e, NULL, NULL, GrB_MAX_FIRST_SEMIRING_FP32, A, t, GrB_DESC_RT0));
  GRB_TRY(GrB_assign(R, NULL, init_flow, e, GrB_ALL, n, S, GrB_DESC_T1));
  GrB_Vector_free(&t); //no longer needed
  

  //begin algorithm loop
  //check if all values of e are zero except the sink index
  GrB_Index n_active;
  GRB_TRY(GrB_Vector_nvals(&n_active, e));

  GrB_Vector J_vector = NULL;
  GRB_TRY(GrB_Vector_new(&J_vector, GrB_INT32, n));

  GrB_Vector y;
  GRB_TRY(GrB_Vector_new(&y, GrB_MFResult, n));

  GrB_Vector delta_vector;
  GRB_TRY(GrB_Vector_new(&delta_vector, GrB_FP32, n));

  GrB_Matrix delta_matrix;
  GRB_TRY(GrB_Matrix_new(&delta_matrix, GrB_FP32, n, n));

  GrB_BinaryOp GrB_MFAddResidual;
  GRB_TRY(GxB_BinaryOp_new(&GrB_MFAddResidual, F_BINARY(MF_updateFlow), GrB_ResidualEdge, GrB_ResidualEdge, GrB_FP32, "MF_updateFlow", GRB_UPDATE_STR));

  GrB_Vector flow;
  GRB_TRY(GrB_Vector_new(&flow, GrB_FP32, n));

  GrB_BinaryOp GrB_MFUpdateHeight;
  GRB_TRY(GxB_BinaryOp_new(&GrB_MFUpdateHeight, F_BINARY(MF_updateHeight), GrB_INT32, GrB_INT32, GrB_MFResult, "MF_updateHeight", MF_UPDATE_HEIGHT_STR));
  
  while(n_active > 0){
    //y<e, struct> = Rxd
    //TODO: modify the index binary op for relabeling and move allocation to the loop.
    GRB_TRY(GrB_mxv(y, e, NULL, Rxd_semiring, R, d, GrB_DESC_RS));
    //update d
    GRB_TRY(GrB_assign(d, NULL, GrB_MFUpdateHeight, y, GrB_ALL, n, NULL));

    //compute delta vector: delta_vec = min(y.r, e)
    GRB_TRY(GrB_eWiseMult(delta_vector, NULL, NULL, GrB_MFeWiseMin, y, e, GrB_DESC_R));
    
    //J_vec = y.j
    GRB_TRY(GrB_apply(J_vector, NULL, NULL, GrB_prune, y, GrB_DESC_R));

    int * jvec_raw;
    float * delta_raw;
    GrB_Index * I;
    GrB_Index * J;
    //extract tuples
    GRB_TRY(GrB_Vector_extractTuples(I, delta_raw, &n, delta_vector));
    GRB_TRY(GrB_Vector_extractTuples(J, jvec_raw, &n, J_vector));

    //build delta matrix
    GRB_TRY(GrB_Matrix_build(delta_matrix, I, (void*) jvec_raw, delta_raw, n, GxB_IGNORE_DUP));
    // deltas = delta - delta^T
    GRB_TRY(GxB_eWiseUnion(delta_matrix, NULL, NULL, GrB_MINUS_FP32, delta_matrix, 0, delta_matrix, 0, GrB_DESC_T1));

    //TODO: update R with eWiseUnion op
    GRB_TRY(GxB_eWiseUnion(R, NULL, NULL, GrB_MFAddResidual, R, 0, delta_matrix, 0, NULL));

    //TODO: update n_active
    GRB_TRY(GrB_reduce(flow, NULL, NULL, GrB_PLUS_FP32, delta_matrix, GrB_DESC_R));
    GRB_TRY(GxB_eWiseUnion(e, NULL, NULL, GrB_PLUS_FP32, e, 0, flow, 0, NULL));
    
    //clear delta matrix
    GRB_TRY(GrB_Matrix_clear(delta_matrix));

    //TODO: free semiring and index binary op
    
  }
  
  LG_FREE_ALL;
  return 0;
}
