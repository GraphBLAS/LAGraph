//************************HOUSE RULES *******************
// user defined typedefs are MF for Max Flow
// when registered with GraphBLAS, replace MF with GrB
// macros are all caps
// op params will be z, y, x accordingly

#include <LAGraphX.h>
#include "LG_internal.h"
#include <LAGraph.h>

#undef LG_FREE_WORK
#undef LG_FREE_ALL

#define LG_FREE_WORK                    \
{                                      \
  GrB_free(&GrB_FlowEdge);             \
  GrB_free(&GrB_CompareTuple);         \
  GrB_free(&GrB_ResultTuple);          \
  GrB_free(&e); \
  GrB_free(&d);\
  GrB_free(&theta);\
  GrB_free(&R);\
  GrB_free(&delta);\
  GrB_free(&e_dup);\
  GrB_free(&A);	   \
  GrB_free(&d_dup);\
  GrB_free(&delta);\
  GrB_free(&delta_vec); \
  GrB_free(&delta_mat); \
  GrB_free(&R_temp1);\
  GrB_free(&R_temp2);\
  GrB_free(&e_temp);\
  GrB_free(&active_set);\
  GrB_free(&map);\
  GrB_free(&y);\
  GrB_free(&yd);\
  GrB_free(&mask_vector);\
  GrB_free(&Jvec);\
  GrB_free(&R_dup);\
  GrB_free(&e_dup);\
  GrB_free(&GrB_UpdateFlows);\
  GrB_free(&GrB_UpdateHeight);\
  GrB_free(&GrB_extractFlows);\
  GrB_free(&GrB_MxeIndexMult);\
  GrB_free(&GrB_MxeMult);\
  GrB_free(&GrB_MxeAdd);\
  GrB_free(&GrB_MxeAddMonoid);\
  GrB_free(&GrB_MxeSemiring);\
  GrB_free(&GrB_extractJ);\
  GrB_free(&GrB_CreateCompareVec);\
  GrB_free(&GrB_RxdSemiring);\
  GrB_free(&GrB_RxdAdd);\
  GrB_free(&GrB_RxdAddMonoid);\
  GrB_free(&GrB_RxdIndexMult);\
  GrB_free(&GrB_RxdMult);\
  GrB_free(&GrB_InitFlows);\
  GrB_free(&GrB_CreateResidual);\
}


#define LG_FREE_ALL \
{ \
  LG_FREE_WORK; \
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
#define GRB_FLOWEDGE_STR "typedef struct{ \nfloat flow; \nfloat capacity; \n} MF_flowEdge;"
#define GRB_RESULTTUPLE_STR                                                    \
  "typedef struct{\nfloat residual;\nint d;\n GrB_Index j;\n} "            \
  "MF_resultTuple;"
#define GRB_COMPARETUPLE_STR "typedef struct{\nfloat residual;\nint di; \nint y_dmin; \nGrB_Index j;\n} MF_compareTuple;"
#define GRB_CR_STR "void MF_CreateResidual(MF_flowEdge *z, const float *y) {\nz->flow = 0;\nz->capacity = (*y);\n}"
#define GRB_RXDMULT_STR "void MF_RxdMult(MF_resultTuple *z, const MF_flowEdge *y, GrB_Index ix, GrB_Index jx, const int *x, GrB_Index iy, GrB_Index jy, const int* theta) {\nfloat r = y->capacity - y->flow;\nif(r > 0){\nz->residual = r;\nz->d = *x;\nz->j = jx;\n}\nelse{\nz->residual = 0;\nz->d = INT32_MAX;\nz->j = -1;\n}\n}"
#define GRB_RXDADD_STR "void MF_RxdAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x) {\nif(y->d < x->d){\nmemcpy(z, y, sizeof(MF_resultTuple));\n}\nelse if(y->d > x->d){\nmemcpy(z, y, sizeof(MF_resultTuple));\n}\nelse{\nif(y->residual > x->residual){\nmemcpy(z, y, sizeof(MF_resultTuple));\n}\nelse if(y->residual < x->residual){\nmemcpy(z, x, sizeof(MF_resultTuple));\n}\nelse{\nif(y->j > x->j){\nmemcpy(z, y, sizeof(MF_resultTuple));\n}\nelse{\nmemcpy(z, x, sizeof(MF_resultTuple));\n}\n}\n}\n}"
#define GRB_INIT_FLOW_STR "void MF_initFlows(MF_flowEdge * z, const MF_flowEdge * y, const float * x){\nz->flow = *x;\n}"
#define GRB_CREATECOMPVEC_STR "void MF_CreateCompareVec(MF_compareTuple *z, const MF_resultTuple *y, const int *x) {\nz->di = (*x);\nz->j = y->j;\nz->residual = y->residual;\nz->y_dmin = y->d;\n}"
#define GRB_EXTRACTJ_STR "void MF_extractJ(int *z, const MF_compareTuple *y) {\nif(y->di < INT32_MAX){\n(*z) = y->j;\n}\n}"
#define GRB_MXEMULT_STR "void MF_MxeMult(MF_resultTuple * z, const MF_compareTuple * y, GrB_Index iy, GrB_Index jy, const float * x, GrB_Index ix, GrB_Index jx, const int* theta){\nif((*x) == 0){\nif(y->di < y->y_dmin-1 || y->di == y->y_dmin+1){\nz->d = y->y_dmin;\nz->residual = y->residual;\nz->j = y->j;\n}\nelse if(y->di == y->y_dmin){ \nif(iy < ix){\nz->d = y->y_dmin;\nz->residual = y->residual;\nz->j = y->j;\n}\n}\n}\nelse{\nz->d = y->y_dmin;\nz->residual = y->residual;\nz->j = y->j;\n}\n}"
#define GRB_MXEADD_STR "void MF_MxeAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x){\nif(x){\nmemcpy(z, x, sizeof(MF_resultTuple));\n}\nelse{\nmemcpy(z, y, sizeof(MF_resultTuple));\n}\n}"
#define GRB_EXTRACTFLOW_STR "void MF_extractFlow(float *z, const MF_resultTuple *y) {\n(*z) = y->residual;\n}"
#define GRB_UPDATEHEIGHT_STR "void MF_updateHeight(int *z, const int *y, const MF_resultTuple *x) {\nif((*y) != x->d+1){\n(*z) = x->d + 1;\n}\nelse{\n(*z) = x->d;\n}\n}"
#define GRB_UPDATEFLOWS_STR "void MF_updateFlow(MF_flowEdge *z, const MF_flowEdge *y, const float *x) {\nz->flow = y->flow + (*x);\n}"

//custom types
typedef struct{
  float flow;
  float capacity;
} MF_flowEdge;

typedef struct{
  float residual;
  int d;
  GrB_Index j;
} MF_resultTuple;

typedef struct{
  float residual;
  int di;
  int y_dmin;
  GrB_Index j;
} MF_compareTuple;
 
void MF_CreateResidual(MF_flowEdge *z, const float *y) {
  z->flow = 0;
  z->capacity = (*y);
}

void MF_RxdMult(MF_resultTuple *z, const MF_flowEdge *y, GrB_Index ix, GrB_Index jx, const int *x, GrB_Index iy, GrB_Index jy, const int* theta) {
  float r = y->capacity - y->flow;
  if(r > 0){
    z->residual = r;
    z->d = *x;
    z->j = jx;
  }
  else{
    z->residual = 0;
    z->d = INT32_MAX;
    z->j = -1;
  }
}

void MF_RxdAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x) {
  if(y->d < x->d){
    memcpy(z, y, sizeof(MF_resultTuple));
  }
  else if(y->d > x->d){
    memcpy(z, y, sizeof(MF_resultTuple));
  }
  else{
    if(y->residual > x->residual){
      memcpy(z, y, sizeof(MF_resultTuple));
    }
    else if(y->residual < x->residual){
      memcpy(z, x, sizeof(MF_resultTuple));
    }
    else{
      if(y->j > x->j){
	memcpy(z, y, sizeof(MF_resultTuple));
      }
      else{
	memcpy(z, x, sizeof(MF_resultTuple));
      }
    }
  }
}


void MF_extractFlow(float *z, const MF_resultTuple *y) {
  (*z) = y->residual;
}

void MF_updateFlow(MF_flowEdge *z, const MF_flowEdge *y, const float *x) {
  z->flow = y->flow + (*x);
}

void MF_updateHeight(int *z, const int *y, const MF_resultTuple *x) {
  if((*y) != x->d+1){
    (*z) = x->d + 1;
  }
  else{
    (*z) = x->d;
  }
}


void MF_extractJ(int *z, const MF_compareTuple *y) {
  if(y->di < INT32_MAX){
    (*z) = y->j;
  }
}

void MF_initFlows(MF_flowEdge * z, const MF_flowEdge * y, const float * x){
  z->flow = *x;
}

void MF_MxeMult(MF_resultTuple * z, const MF_compareTuple * y, GrB_Index iy, GrB_Index jy, const float * x, GrB_Index ix, GrB_Index jx, const int* theta){
  if((*x) == 0){
    if(y->di < y->y_dmin-1 || y->di == y->y_dmin+1){
      z->d = y->y_dmin;
      z->residual = y->residual;
      z->j = y->j;
    }
    else if(y->di == y->y_dmin){ //check this
      if(iy < ix){
	z->d = y->y_dmin;
	z->residual = y->residual;
	z->j = y->j;
      }
    }
  }
  else{
    z->d = y->y_dmin;
    z->residual = y->residual;
    z->j = y->j;
  }
}

void MF_MxeAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x){
  if(x){
    memcpy(z, x, sizeof(MF_resultTuple));
  }
  else{
    memcpy(z, y, sizeof(MF_resultTuple));
  }
}

void MF_CreateCompareVec(MF_compareTuple *z, const MF_resultTuple *y, const int *x) {
  z->di = (*x);
  z->j = y->j;
  z->residual = y->residual;
  z->y_dmin = y->d;
}

int LAGraph_MaxFlow(LAGraph_Graph G, GrB_Index* S, GrB_Index *T, int * f, char *msg){

  //plan of ATTACK ****************************************
  //1. Create R as an anti-symmetric graph of flow edges from Adj of G, using ewise union
  //
  //2. initialize d and e, where d is dense and e contains the saturated flows of the first level of the graph, do select mxv for e
  //
  //
  //3. do a column assign to place the saturated flows in the R matrix at columns index S, do a transpose iun descriptor
  //
  //
  //4. begin execution loop, do y<e, struct> = Rxd
  //5. take y and d to create vector yd of type compareTuple
  //6. unpack the vector and create map matrix
  //7. compute y = map * e
  //8. update d from y
  //9. compute delta from min(y.r, e)
  //10. create antisymmetric delta matrix
  //
  //
  //11. update R
  //12. reduce delta matrix to df vector and add to e
  //13. get number of active nodes from e through an extract op
  //
  // 14. set f to value of e(T)

  //types
  GrB_Type GrB_FlowEdge;
  GrB_Type GrB_ResultTuple;
  GrB_Type GrB_CompareTuple;

  //to create R
  GrB_UnaryOp GrB_CreateResidual;
  GrB_Matrix A = G->A;
  GrB_Index n;
  GrB_Matrix R_temp1, R_temp2, R;

  //to init R with initial saturated flows
  GrB_Vector e_temp;
  GrB_Vector e;
  GrB_BinaryOp GrB_InitFlows;

  //create height vector
  GrB_Vector d;
 

  //active_set and n_active
  GrB_Vector active_set;
  GrB_Vector mask_vector;
  GrB_Index n_active;

  //semiring and vectors for y<e, struct> = R x d
  GrB_Vector y;
  GzB_IndexBinaryOp GrB_RxdIndexMult;
  GrB_BinaryOp GrB_RxdAdd, GrB_RxdMult;
  GrB_Monoid GrB_RxdAddMonoid;
  GrB_Semiring GrB_RxdSemiring;
  GrB_Scalar theta;
 
  //binary op and yd
  GrB_Vector yd;
  GrB_BinaryOp GrB_CreateCompareVec;
  
  //utility vectors, Matrix, and ops for mapping
  GrB_Matrix map;
  GrB_Vector Jvec;
  GrB_UnaryOp GrB_extractJ;
  GrB_Index* Jmap;
  GrB_Index* Imap;
  GrB_Index* Jvec_value;
  MF_compareTuple* yd_value;
 
  //map x e semiring
  GrB_Semiring GrB_MxeSemiring;
  GrB_Monoid GrB_MxeAddMonoid;
  GrB_BinaryOp GrB_MxeAdd, GrB_MxeMult;
  GzB_IndexBinaryOp GrB_MxeIndexMult;

  //residual flow vec
  GrB_Vector residual_vec;
  GrB_UnaryOp GrB_extractFlows;
 
  //delta structures
  GrB_Vector delta_vec;
  GrB_Matrix delta, delta_mat;
  GrB_Index* Idelta;
  float* delta_raw;

  //relabel
  GrB_Vector d_dup;

  //update height
  GrB_BinaryOp GrB_UpdateHeight;

  //update R structure
  GrB_Matrix R_dup;
  GrB_BinaryOp GrB_UpdateFlows;

  //update e
  GrB_Vector e_dup;
  
  //create types for computation
  GRB_TRY(GxB_Type_new(&GrB_FlowEdge, sizeof(MF_flowEdge), "MF_FlowEdge", GRB_FLOWEDGE_STR));
  GRB_TRY(GxB_Type_new(&GrB_ResultTuple, sizeof(MF_resultTuple), "MF_resultTuple", GRB_RESULTTUPLE_STR));
  GRB_TRY(GxB_Type_new(&GrB_CompareTuple, sizeof(MF_compareTuple), "MF_compareTuple", GRB_COMPARETUPLE_STR));

  //create R
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidual, F_UNARY(MF_CreateResidual), GrB_FlowEdge , GrB_FP32, "MF_CreateResidual", GRB_CR_STR));
  GRB_TRY(GrB_Matrix_nrows(&n, A));
  GRB_TRY(GrB_Matrix_new(&R_temp1, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_Matrix_new(&R_temp2, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_Matrix_new(&R, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_apply(R_temp1, NULL, NULL, GrB_CreateResidual, A, NULL));
  GRB_TRY(GrB_apply(R_temp2, NULL, NULL, GrB_CreateResidual, A, GrB_DESC_T0));
  GRB_TRY(GrB_assign(R, NULL, NULL, R_temp1, GrB_ALL, n, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(R, A, NULL, R_temp2, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
  GrB_free(&R_temp2);
  GrB_free(&R_temp1);

  //init R with initial saturated flows
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitFlows, F_BINARY(MF_initFlows), GrB_FlowEdge, GrB_FlowEdge, GrB_FP32, "MF_initFlows", GRB_INIT_FLOW_STR));
  GRB_TRY(GrB_Vector_new(&e_temp, GrB_FP32, n));
  GRB_TRY(GrB_Vector_new(&e, GrB_FP32, n));
  GRB_TRY(GrB_assign(e_temp, NULL, NULL, 1, S, 1, NULL));
  GRB_TRY(GrB_mxv(e, NULL, NULL, GxB_MAX_FIRST_FP32, A, e_temp, NULL));
  GRB_TRY(GrB_assign(R, NULL, GrB_InitFlows, e, GrB_ALL, n, *S, GrB_DESC_T1));//transpose R?
  GrB_free(&e_temp);

  //create and init d vector
  GRB_TRY(GrB_Vector_new(&d, GrB_INT32, n));
  GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(d, NULL, NULL, n, S, 1, NULL));

  //extract n_active from e masking T and S then assign to e
  GRB_TRY(GrB_Vector_new(&active_set, GrB_FP32, n));
  GRB_TRY(GrB_Vector_new(&mask_vector, GrB_BOOL, n)); //keep as bool?
  GRB_TRY(GrB_assign(mask_vector, NULL, NULL, true, T, 1, NULL));
  GRB_TRY(GrB_assign(mask_vector, NULL, NULL, true, S, 1, NULL));
  GRB_TRY(GrB_select(active_set, mask_vector, NULL, GrB_VALUEGE_FP32, e, 0, GrB_DESC_RSC));
  GRB_TRY(GrB_assign(e, NULL, NULL, active_set, GrB_ALL, n, GrB_DESC_SC));
  GRB_TRY(GrB_Vector_nvals(&n_active, active_set));

  //create semiring and vectors for y<e, struct> = R x d
  GRB_TRY(GrB_Scalar_new(&theta, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement_INT32(theta, 0));
  GRB_TRY(GrB_Vector_new(&y, GrB_ResultTuple, n));
  GRB_TRY(GzB_IndexBinaryOp_new(&GrB_RxdIndexMult, F_INDEX_BINARY(MF_RxdMult), GrB_ResultTuple, GrB_FlowEdge, GrB_INT32, GrB_INT32, "MF_RxdMult", GRB_RXDMULT_STR));
  GRB_TRY(GzB_BinaryOp_new_IndexOp(&GrB_RxdMult, GrB_RxdIndexMult, theta));
  GRB_TRY(GxB_BinaryOp_new(&GrB_RxdAdd, F_BINARY(MF_RxdAdd), GrB_ResultTuple, GrB_ResultTuple, GrB_ResultTuple, "MF_RxdAdd", GRB_RXDADD_STR));
  MF_resultTuple id = {.d = INT32_MAX, .j = -1, .residual = 0};
  GRB_TRY(GrB_Monoid_new_UDT(&GrB_RxdAddMonoid, GrB_RxdAdd, &id));
  GRB_TRY(GrB_Semiring_new(&GrB_RxdSemiring, GrB_RxdAddMonoid, GrB_RxdMult));

  //create binary op and yd
  GRB_TRY(GrB_Vector_new(&yd, GrB_CompareTuple, n));
  GRB_TRY(GxB_BinaryOp_new(&GrB_CreateCompareVec, F_BINARY(MF_CreateCompareVec), GrB_CompareTuple, GrB_ResultTuple, GrB_INT32, "MF_CreateCompareVec", GRB_CREATECOMPVEC_STR));

  //create utility vectors, Matrix, and ops for mapping
  GRB_TRY(GrB_Vector_new(&Jvec, GrB_INT32, n));
  GRB_TRY(GrB_Matrix_new(&map, GrB_CompareTuple, n,n));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractJ, F_UNARY(MF_extractJ), GrB_INT32, GrB_CompareTuple, "MF_extractJ", GRB_EXTRACTJ_STR));

  //create map x e semiring
  GRB_TRY(GzB_IndexBinaryOp_new(&GrB_MxeIndexMult, F_INDEX_BINARY(MF_MxeMult), GrB_ResultTuple, GrB_CompareTuple, GrB_FP32, GrB_INT32, "MF_MxeMult", GRB_MXEMULT_STR));
  GRB_TRY(GzB_BinaryOp_new_IndexOp(&GrB_MxeMult, GrB_MxeIndexMult, theta));
  GRB_TRY(GxB_BinaryOp_new(&GrB_MxeAdd, F_BINARY(MF_MxeAdd), GrB_ResultTuple, GrB_ResultTuple, GrB_ResultTuple, "MF_MxeAdd", GRB_MXEADD_STR));
  GRB_TRY(GrB_Monoid_new_UDT(&GrB_MxeAddMonoid, GrB_MxeAdd, &id));
  GRB_TRY(GrB_Semiring_new(&GrB_MxeSemiring, GrB_MxeAddMonoid, GrB_MxeMult));

  //create flow vec
  GRB_TRY(GrB_Vector_new(&residual_vec, GrB_FP32, n));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractFlows, F_UNARY(MF_extractFlow), GrB_FP32, GrB_ResultTuple, "MF_extractFlow", GRB_EXTRACTFLOW_STR));

  GRB_TRY(GrB_Matrix_new(&delta_mat, GrB_FP32, n, n));
  GRB_TRY(GrB_Matrix_new(&delta, GrB_FP32, n, n));
  GRB_TRY(GrB_Vector_new(&delta_vec, GrB_FP32, n));

  //relable structures
  GRB_TRY(GrB_Vector_new(&d_dup, GrB_INT32, n));

  //update height binary op
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateHeight, F_BINARY(MF_updateHeight), GrB_INT32, GrB_INT32, GrB_ResultTuple, "MF_updateHeight", GRB_UPDATEHEIGHT_STR));

  //update R structure
  GRB_TRY(GrB_Matrix_new(&R_dup, GrB_FlowEdge, n, n));
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateFlows, F_BINARY(MF_updateFlow), GrB_FlowEdge, GrB_FlowEdge, GrB_FP32, "MF_updateflows", GRB_UPDATEFLOWS_STR));

  //update e structures
  GRB_TRY(GrB_Vector_new(&e_dup, GrB_FP32, n));
  
  while(n_active > 0){
    
    //y<e, struct> = R x d
    GRB_TRY(GrB_mxv(y, e, NULL, GrB_RxdSemiring, R, d, GrB_DESC_RS));

    //create yd vector of type compare tuple
    GRB_TRY(GrB_eWiseMult(yd, NULL, NULL, GrB_CreateCompareVec, y, d, GrB_DESC_R));

    //create map matrix from yd
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractJ, yd, GrB_DESC_R));
    GRB_TRY(GrB_Vector_extractTuples(Jmap, Jvec_value, &n, Jvec));
    GRB_TRY(GrB_Vector_extractTuples(Imap, (void*)yd_value, &n, yd));
    GRB_TRY(GrB_Matrix_build(map, Imap, Jvec_value, (void*)yd_value, n, GxB_IGNORE_DUP));

    //y = map x e
    GRB_TRY(GrB_mxv(y, NULL, NULL, GrB_MxeSemiring, map, e, GrB_DESC_R));

    //relable, update heights
    GRB_TRY(GrB_Vector_dup(&d_dup, d));
    GRB_TRY(GxB_eWiseUnion(d, NULL, NULL, GrB_UpdateHeight, d_dup, NULL, y, NULL, NULL));

    //extract residual flows from y
    GRB_TRY(GrB_apply(residual_vec, NULL, NULL, GrB_extractFlows, y, GrB_DESC_R));

    //.min(flow_vec and e)
    GRB_TRY(GrB_eWiseMult(delta_vec, NULL, NULL, GrB_MIN_FP32, residual_vec, e, GrB_DESC_R));

    //create delta matrix from delta vector
    GRB_TRY(GrB_Vector_extractTuples(Idelta, delta_raw, &n, delta_vec));
    GRB_TRY(GrB_Matrix_build(delta, Idelta, Jvec_value, delta_raw, n, GxB_IGNORE_DUP));

    //make delta anti-symmetric
    GRB_TRY(GxB_eWiseUnion(delta_mat, NULL, NULL, GrB_MINUS_FP32, delta, NULL, delta, NULL, GrB_DESC_RT1));

    //update R
    GRB_TRY(GrB_Matrix_dup(&R_dup, R));
    GRB_TRY(GrB_eWiseMult(R, NULL, NULL, GrB_UpdateFlows, R_dup, delta_mat, NULL));

    //reduce R to delta_vec
    GRB_TRY(GrB_reduce(delta_vec, NULL, NULL, GrB_PLUS_FP32, delta_mat, GrB_DESC_R));

    //add to e
    GRB_TRY(GrB_Vector_dup(&e_dup, e));
    GRB_TRY(GxB_eWiseUnion(e, NULL, NULL, GrB_PLUS_FP32, e_dup, NULL, delta_vec, NULL, NULL));

    //extract active nodes
    GRB_TRY(GrB_select(active_set, mask_vector, NULL, GrB_VALUEGE_FP32, e, 0, GrB_DESC_RSC));
    GRB_TRY(GrB_assign(e, NULL, NULL, active_set, GrB_ALL, n, GrB_DESC_SC));
    GRB_TRY(GrB_Vector_nvals(&n_active, active_set));

    //clear map and delta
    GRB_TRY(GrB_Matrix_clear(map));
    GRB_TRY(GrB_Matrix_clear(delta));
    
  }

  //set f
  GRB_TRY(GrB_Vector_extractElement(f, e, *T));
  
  LG_FREE_WORK;
  return GrB_SUCCESS;
}
