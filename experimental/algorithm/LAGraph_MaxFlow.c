
#include <LG_internal.h>
#include <LAGraphX.h>
#include <LAGraph.h>

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

// casting for binary op
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

// strings for JIT
#define GRB_FLOWEDGE_STR "Typedef struct{ float flow; float capacity; } GrB_Flow_Edge;"
#define GRB_MFRESULT_STR "typedef struct{ float residual; int d; GrB_Index j; } MF_result_tuple; "
#define GRB_CRF_STR "void CreateResidualForward_UOp(GrB_Flow_Edge *f, const float *cap) {f->flow = 0; f->capacity = (*cap);}"
#define GRB_CRB_STR "void CreateResidualBackward_UOp(GrB_Flow_Edge *f, const float *cap) {f->flow = 0; f->capacity = (*cap);}"
#define GZB_MULT_STR "void Rxd_MultBOp(MF_result_tuple *y, const GrB_Flow_Edge *R, GrB_Index ix, GrB_Index jx, const int *d, GrB_Index iy, GrB_Index ik, const int* theta) { float r = R->capacity - R->flow; if(r > 0){ y->residual = r; y->d = *d; y->j = jx; } else{ y->residual = 0; y->d = INT32_MAX; y->j = -1; } }"
#define GRB_ADD_STR "void Rxd_AddMonoid(MF_result_tuple * z, const MF_result_tuple * x, const MF_result_tuple * y) {if(x->d < y->d){ memcpy(z, x, sizeof(MF_result_tuple)); } else if(x->d > y->d){ memcpy(z, y, sizeof(MF_result_tuple));}else{if(x->residual > y->residual){ memcpy(z, x, sizeof(MF_result_tuple)); } else if(x->residual < y->residual){ memcpy(z, y, sizeof(MF_result_tuple)); } else{ if(x->j > y->j){ memcpy(z, x, sizeof(MF_result_tuple)); } else{ memcpy(z, y, sizeof(MF_result_tuple)); } } } }"

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
void CreateResidualForward_UOp(GrB_Flow_Edge *f, const float *cap) {
  f->flow = 0;
  f->capacity = (*cap);
}

//Make this an Unary apply op
//takes capacity from adj and creates flows for the back edges 
void CreateResidualBackward_UOp(GrB_Flow_Edge *f, const float *cap) {
  f->flow = 0;
  f->capacity = (*cap); //change to INF ??
}

//TO-DO: create add operation for flow edges.

void Rxd_MultBOp(MF_result_tuple *y, const GrB_Flow_Edge *R, GrB_Index ix, GrB_Index jx, const int *d, GrB_Index iy, GrB_Index ik, const int* theta) {
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

// R is resulting residual graph
// f is max flow
int LAGraph_MaxFlow(LAGraph_Graph G, GrB_Index S, GrB_Index T, int * f, char *msg){
  
  //create semirings and types
  GrB_Type GrB_ResidualEdge;
  GRB_TRY(GxB_Type_new(&GrB_ResidualEdge, sizeof(GrB_Flow_Edge), "GrB_Flow_Edge", GRB_FLOWEDGE_STR));

  GrB_Type GrB_MFResult;
  GRB_TRY(GxB_Type_new(&GrB_MFResult, sizeof(MF_result_tuple), "MF_result_tuple", GRB_MFRESULT_STR));

  GrB_UnaryOp GrB_CRF_UOp;
  GRB_TRY(GxB_UnaryOp_new(&GrB_CRF_UOp, F_UNARY(CreateResidualForward_UOp), GrB_ResidualEdge, GrB_FP32, "CreateResidualForward_UOp", GRB_CRF_STR));

  GrB_UnaryOp GrB_CRB_UOp;
  GRB_TRY(GxB_UnaryOp_new(&GrB_CRB_UOp, F_UNARY(CreateResidualBackward_UOp), GrB_ResidualEdge, GrB_FP32, "CreateResidualForward_UOp", GRB_CRB_STR));

  GzB_IndexBinaryOp GzB_MFMult;
  GrB_BinaryOp GrB_RxdMult;
  GrB_Scalar unused_theta;
  GRB_TRY(GrB_Scalar_new(&unused_theta, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement_INT32(unused_theta, 0));
  GRB_TRY(GzB_IndexBinaryOp_new(&GzB_MFMult, F_INDEX_BINARY(Rxd_MultBOp), GrB_MFResult, GrB_ResidualEdge, GrB_INT32, GrB_INT32,  "Rxd_MultBOp", GZB_MULT_STR));
  GRB_TRY(GzB_BinaryOp_new_IndexOp(&GrB_RxdMult, GzB_MFMult, unused_theta));

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
  GRB_TRY(GrB_apply(R, NULL, NULL, GrB_CRF_UOp, A, NULL));
  GRB_TRY(GrB_apply(R, NULL, NULL, GrB_CRB_UOp, A, GrB_DESC_T1));
  
  //create d (height) vector and e (excess) vector
  GrB_Vector d = NULL;
  GrB_Vector e = NULL;
  GRB_TRY(GrB_Vector_new(&e, GrB_FP32, n));
  GRB_TRY(GrB_Vector_new(&d, GrB_INT32, n));

  //init e and d
  GrB_Scalar size;
  GRB_TRY(GrB_Scalar_new(&size, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement_INT32(size, n));
  GRB_TRY(GrB_Vector_setElement(d, size, S)); 
  
  LG_FREE_ALL;
  return 0;
}
