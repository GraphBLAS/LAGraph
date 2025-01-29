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
  GrB_free(&active_set);\
  GrB_free(&map);\
  GrB_free(&y);\
  GrB_free(&yd);\
  GrB_free(&mask_vector);\
  GrB_free(&Jvec);\
  GrB_free(&GrB_Prune);\
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
  GrB_free(&GrB_InitForwardFlows);\
  GrB_free(&GrB_InitBackwardFlows);\
  GrB_free(&GrB_CreateResidualForward);\
  GrB_free(&GrB_CreateResidualBackward);\
  GrB_free(&zero_int32);\
  GrB_free(&zero_fp32);\
  GrB_free(&Re);\
  GrB_free(&invariant);\
  GrB_free(&GrB_InvariantCheck);\
  GrB_free(&check);\
  GrB_free(&GrB_extractYJ);\
}

#define LG_FREE_ALL \
{ \
  LG_FREE_WORK; \
}

#define LEN 512

//casting for unary ops
#define F_UNARY(f) ((void (*)(void *, const void *))f)

// casting for index binary ops
#define F_INDEX_BINARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void *, GrB_Index, GrB_Index, const void *)) f)

// casting for index unary op
//#define F_INDEX_UNARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void*)) f)

// casting for binary op
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

// strings for JIT
#define GRB_EXTRACTYJ_STR "void MF_extractYJ(int *z, const MF_resultTuple *y) {" \
  "(*z) = y->j;" \
"}"

#define GRB_GETRES_STR "void MF_getResidual(double * z, const MF_flowEdge * y){" \
  "*z = y->capacity - y->flow;" \
"}"

#define GRB_PRUNE_STR "void MF_Prune(bool * z, const MF_resultTuple * y, GrB_Index iy, GrB_Index jy, const int * theta){"\
  "if(y->j != *theta){"\
    "*z = true;" \
  "}" \
  "else{" \
    "*z = false;" \
  "}" \
"}"

#define GRB_FLOWEDGE_STR "typedef struct{"\
  "double flow;"\
  "double capacity;"\
"} MF_flowEdge;"

#define GRB_RESULTTUPLE_STR "typedef struct{"\
  "double residual;"\
  "int d;"\
  "GrB_Index j;"\
"} MF_resultTuple;"

#define GRB_COMPARETUPLE_STR                                                   \
  "typedef struct{"                                                            \
  "double residual;"                                                            \
  "int di;"                                                                    \
  "int y_dmin;"                                                                \
  "GrB_Index j;"                                                               \
  "} MF_compareTuple;"

#define GRB_CRF_STR "void MF_CreateResidualForward(MF_flowEdge *z, const double *y) {"\
  "z->flow = 0;"\
  "z->capacity = (*y);"\
"}"

#define GRB_CRB_STR "void MF_CreateResidualBackward(MF_flowEdge *z, const double *y)"\
 "{"\
  "z->flow " \
  "= 0;\nz->capacity = 0;"\
"}"

#define GRB_RXDMULT_STR "void MF_RxdMult(MF_resultTuple *z, const MF_flowEdge *y, GrB_Index iy, GrB_Index jy, const int *x, GrB_Index ix, GrB_Index jx, const int* theta) {"\
  "double r = y->capacity - y->flow;"\
  "if(r > 0){"\
    "z->residual = r;"\
    "z->d = *x;"\
    "z->j = jy;"\
  "}"\
  "else{"\
    "z->residual = 0;"\
    "z->d = INT32_MAX;"\
    "z->j = -1;"\
  "}"\
"}"


#define GRB_RXDADD_STR "void MF_RxdAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x) {"\
  "if(y->d < x->d){"\
    "memcpy(z, y, sizeof(MF_resultTuple));"\
  "}"\
  "else if(y->d > x->d){"\
    "memcpy(z, x, sizeof(MF_resultTuple));"\
  "}"\
  "else{"\
    "if(y->residual > x->residual){"\
      "memcpy(z, y, sizeof(MF_resultTuple));"\
    "}"\
    "else if(y->residual < x->residual){"\
      "memcpy(z, x, sizeof(MF_resultTuple));"\
    "}"\
    "else{"\
      "if(y->j > x->j){"\
	"memcpy(z, y, sizeof(MF_resultTuple));"\
      "}"\
      "else{"\
	"memcpy(z, x, sizeof(MF_resultTuple));"\
      "}"\
    "}"\
  "}"\
"}" 

#define GRB_INITFLOWF_STR "void MF_initFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){"\
  "z->flow = x->flow + y->flow;"\
  "z->capacity = y->capacity;"\
"}"

#define GRB_INITFLOWB_STR "void MF_initFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){"\
  "z->flow = y->flow - x->flow;"\
  "z->capacity = y->capaciy;"\
"}"

#define GRB_CREATECOMPVEC_STR "void MF_CreateCompareVec(MF_compareTuple *z, const MF_resultTuple *y, const int *x) {"\
  "z->di = (*x);"\
  "z->j = y->j;"\
  "z->residual = y->residual;"\
  "z->y_dmin = y->d;"\
"}"

#define GRB_EXTRACTJ_STR "void MF_extractJ(int *z, const MF_compareTuple *y) {"\
  "if(y->j != -1){"\
    "(*z) = y->j;"\
  "}"\
"}"

#define GRB_MXEMULT_STR "void MF_MxeMult(MF_resultTuple * z, const MF_compareTuple * y, GrB_Index iy, GrB_Index jy, const double * x, GrB_Index ix, GrB_Index jx, const int* theta){" \
  "if(y->di == y->y_dmin && (*x) > 0){" \
    "if(iy < jy){" \
      "z->d = y->y_dmin;" \
      "z->residual = y->residual;" \
      "z->j = y->j;" \
    "}" \
    "else{"  \
      "z->d = INT32_MAX;" \
      "z->j = -1;" \
      "z->residual = 0;" \
    "}" \
  "}" \
  "else if(y->di == y->y_dmin - 1 && (*x) > 0){" \
    "z->d = INT32_MAX;" \
    "z->j = -1;" \
    "z->residual = 0;" \
  "}" \
  "else if(y->di < y->y_dmin-1 || y->di == y->y_dmin+1 || y->di == y->y_dmin-1 || y->di == y->y_dmin){" \
    "z->d = y->y_dmin;" \
    "z->residual = y->residual;" \
    "z->j = y->j;" \
  "}" \
  "else{" \
    "z->d = INT32_MAX;" \
    "z->residual = 0;" \
    "z->j = -1;" \
  "}" \
"}"


#define GRB_MXEADD_STR "void MF_MxeAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x){"\
  "if(x != NULL){"\
    "memcpy(z, x, sizeof(MF_resultTuple));"\
  "}"\
  "else{"\
    "memcpy(z, y, sizeof(MF_resultTuple));"\
  "}"\
"}"

#define GRB_EXTRACTFLOW_STR "void MF_extractFlow(double *z, const MF_resultTuple *y) {\n(*z) = y->residual;\n}"
#define GRB_UPDATEHEIGHT_STR "void MF_updateHeight(int *z, const int *y, const MF_resultTuple *x) {"\
  "if((*y) < x->d+1){"\
    "(*z) = x->d + 1;"\
  "}"\
  "else if ((*y) == x->d+1){"\
    "(*z) = (*y);"\
  "}"\
"}"

#define GRB_UPDATEFLOWS_STR "void MF_updateFlow(MF_flowEdge *z, const MF_flowEdge *y, const double *x) {"\
  "z->capacity = y->capacity;"\
  "z->flow = y->flow + (*x);"\
"}"

#define GRB_MAKEF_STR "void MF_MakeFlow(MF_flowEdge * z, const double * y){"\
  "z->capacity = 0;"\
  "z->flow = (*y);"\
"}"

#define GRB_INV_STR "void MF_CheckInvariant(bool *z, const int *y, const MF_resultTuple *x) {"\
  "(*z) = ((*y) == x->d+1);"\
"}"


//custom types
typedef struct{
  double flow;
  double capacity;
} MF_flowEdge;

typedef struct{
  double residual;
  int d;
  GrB_Index j;
} MF_resultTuple;

typedef struct{
  double residual;
  int di;
  int y_dmin;
  GrB_Index j;
} MF_compareTuple;
 
void MF_CreateResidualForward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = (*y);
}

void MF_CreateResidualBackward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = 0;
}

void MF_RxdMult(MF_resultTuple *z, const MF_flowEdge *y, GrB_Index iy, GrB_Index jy, const int *x, GrB_Index ix, GrB_Index jx, const int* theta) {
  double r = y->capacity - y->flow;
  if(r > 0){
    z->residual = r;
    z->d = *x;
    z->j = jy;
  }
  else if(r==0){
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
    memcpy(z, x, sizeof(MF_resultTuple));
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


void MF_extractFlow(double *z, const MF_resultTuple *y) {
  (*z) = y->residual;
}

void MF_updateFlow(MF_flowEdge *z, const MF_flowEdge *y, const double *x) {
  z->capacity = y->capacity;
  z->flow = y->flow + (*x);
}

void MF_updateHeight(int *z, const int *y, const MF_resultTuple *x) {
  if((*y) < x->d+1){
    (*z) = x->d + 1;
  }
  else if ((*y) == x->d+1){
    (*z) = (*y);
  }
}


void MF_extractJ(int *z, const MF_compareTuple *y) {
  (*z) = y->j;
}

void MF_extractYJ(int *z, const MF_resultTuple *y) {
  (*z) = y->j;
}


void MF_initForwardFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){
  z->flow = x->flow + y->flow;
  z->capacity = y->capacity;
}

void MF_initBackwardFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){
  z->flow = y->flow - x->flow;
  z->capacity = y->capacity; 
}

void MF_MxeMult(MF_resultTuple * z, const MF_compareTuple * y, GrB_Index iy, GrB_Index jy, const double * x, GrB_Index ix, GrB_Index jx, const int* theta){
  if(y->di == y->y_dmin && (*x) > 0){ //check this
    if(iy < jy){
      z->d = y->y_dmin;
      z->residual = y->residual;
      z->j = y->j;
    } //add else to populate with empty tuple, prune after.
    else{
      z->d = INT32_MAX;
      z->j = -1;
      z->residual = 0;
    }
  }
  else if(y->di == y->y_dmin - 1 && (*x) > 0){
    z->d = INT32_MAX;
    z->j = -1;
    z->residual = 0;
  }
  else if(y->di < y->y_dmin-1 || y->di == y->y_dmin+1 || y->di == y->y_dmin-1 || y->di == y->y_dmin){
    z->d = y->y_dmin;
    z->residual = y->residual;
    z->j = y->j;
  }
  else{ //change later to signify the removal of the node from the active set since flow cannot be pushed anywhere.
    z->d = INT32_MAX;
    z->residual = 0;
    z->j = -1;
  }
}

void MF_MxeAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x){
  if(x != NULL){
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

void MF_Prune(bool * z, const MF_resultTuple * y, GrB_Index iy, GrB_Index jy, const int * theta){
  if(y->j != *theta){
    *z = true;
  }
  else{
    *z = false;
  }
}

void MF_MakeFlow(MF_flowEdge * z, const double * y){
  z->capacity = 0;
  z->flow = (*y);
}

void print_flowMtx(const GrB_Matrix mtx) {
  GxB_Iterator iter;
  GxB_Iterator_new(&iter);
  GrB_Info info = GxB_Matrix_Iterator_attach(iter, mtx, NULL);
  if(info < 0){
    printf("error with matrix passed in");
  }
  info = GxB_Matrix_Iterator_seek(iter, 0);
  while(info != GxB_EXHAUSTED){
    GrB_Index i, j;
    GxB_Matrix_Iterator_getIndex(iter, &i, &j);
    MF_flowEdge e;
    GxB_Iterator_get_UDT(iter, &e);
    printf("(%ld, %ld)         (capacity :%f, flow: %f) \n", i, j, e.capacity, e.flow);
    info = GxB_Matrix_Iterator_next(iter);
  }
  GrB_free(&iter);
}

void print_MapMtx(const GrB_Matrix mtx) {
  GxB_Iterator iter;
  GxB_Iterator_new(&iter);
  GrB_Info info = GxB_Matrix_Iterator_attach(iter, mtx, NULL);
  if(info < 0){
    printf("error with matrix passed in");
  }
  info = GxB_Matrix_Iterator_seek(iter, 0);
  while(info != GxB_EXHAUSTED){
    GrB_Index i, j;
    GxB_Matrix_Iterator_getIndex(iter, &i, &j);
    MF_compareTuple e;
    GxB_Iterator_get_UDT(iter, &e);
    printf("(%ld, %ld)         (height: %d, y.dmin: %d, J: %ld, residual: %lf) \n", i, j, e.di, e.y_dmin, e.j, e.residual);
    info = GxB_Matrix_Iterator_next(iter);
  }
  GrB_free(&iter);
}

void print_resultVec(const GrB_Vector vec) {
  GxB_Iterator iter;
  GxB_Iterator_new(&iter);
  GrB_Info info = GxB_Vector_Iterator_attach(iter, vec, NULL);
  if(info < 0){
    printf("error with matrix passed in");
  }
  info = GxB_Vector_Iterator_seek(iter, 0);
  while(info != GxB_EXHAUSTED){
    GrB_Index i;
    i = GxB_Vector_Iterator_getIndex(iter);
    MF_resultTuple e;
    GxB_Iterator_get_UDT(iter, &e);
    printf("(%ld, 0)         (height: %d, J: %ld, residual: %f) \n", i, e.d, e.j, e.residual);
    info = GxB_Vector_Iterator_next(iter);
  }
  GrB_free(&iter);
}

void print_compareVec(const GrB_Vector vec) {
  GxB_Iterator iter;
  GxB_Iterator_new(&iter);
  GrB_Info info = GxB_Vector_Iterator_attach(iter, vec, NULL);
  if(info < 0){
    printf("error with matrix passed in");
  }
  info = GxB_Vector_Iterator_seek(iter, 0);
  while(info != GxB_EXHAUSTED){
    GrB_Index i;
    i = GxB_Vector_Iterator_getIndex(iter);
    MF_compareTuple e;
    GxB_Iterator_get_UDT(iter, &e);
    printf("(%ld, 0)         (height: %d, y.dmin: %d, J: %ld, residual: %lf) \n", i, e.di, e.y_dmin, e.j, e.residual);
    info = GxB_Vector_Iterator_next(iter);
  }
  GrB_free(&iter);
}


void MF_CheckInvariant(bool *z, const int *y, const MF_resultTuple *x) {
  (*z) = ((*y) == x->d+1);
}

void MF_getResidual(double * z, const MF_flowEdge * y){
  *z = y->capacity - y->flow;
}

#define GLOBAL_RELABEL                                                         \
  {                                                                            \
    GrB_Vector parent, lvl;                                                    \
    GrB_UnaryOp GrB_GetResidual;                                               \
    GrB_Matrix res_mat, modified_res_mat;                                      \
    LAGraph_Graph res_graph;                                                   \
    GrB_Vector_new(&parent, GrB_INT64, n);				\
    GrB_Vector_new(&lvl, GrB_INT64, n);					\
    GrB_Matrix_new(&res_mat, GrB_FP64, n, n);				\
    GrB_Matrix_new(&modified_res_mat, GrB_FP64, n, n);			\
    GxB_UnaryOp_new(&GrB_GetResidual, F_UNARY(MF_getResidual), GrB_FP64,       \
                    GrB_FlowEdge, "MF_getResidual", GRB_GETRES_STR);           \
    GrB_apply(res_mat, NULL, NULL, GrB_GetResidual, R, NULL);                  \
    GrB_Matrix_dup(&modified_res_mat, res_mat);                                \
    GrB_select(modified_res_mat, NULL, NULL, GrB_VALUEGT_FP64, res_mat, 0,     \
               GrB_DESC_R);                                                          \
    LAGraph_New(&res_graph, &modified_res_mat, LAGraph_ADJACENCY_DIRECTED,  \
                   msg);                                                      \
    LAGraph_Cached_AT(res_graph, msg);                                     \
    LAGraph_Cached_OutDegree(res_graph, msg);                              \
    LAGr_BreadthFirstSearch(&lvl, &parent, res_graph, T, msg);             \
    GrB_assign(d, NULL, NULL, lvl, GrB_ALL, n, GrB_DESC_R);                    \
    GrB_assign(d, lvl, NULL, 0, GrB_ALL, n, GrB_DESC_SC);                      \
    GrB_free(&parent);                                                         \
    GrB_free(&lvl);                                                            \
    GrB_free(&GrB_GetResidual);                                                \
    GrB_free(&res_mat);                                                        \
    GrB_free(&modified_res_mat);                                               \
    LAGraph_Delete(&res_graph, msg);                                       \
  }

//GrB_assign(d, lvl, NULL, 0, GrB_ALL, n, GrB_DESC_SC);                      \
  
int LAGraph_MaxFlow(LAGraph_Graph G, GrB_Index S, GrB_Index T, double * f, char *msg){

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
  //7. compute y = map * e, make sure e is a dense vector with explicit zeros
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
  GrB_UnaryOp GrB_CreateResidualForward, GrB_CreateResidualBackward;
  GrB_Matrix A = G->A;
  GrB_Index n;
  GrB_Matrix R_temp1, R_temp2, R;
  GrB_Matrix_nrows(&n, A);

  //to init R with initial saturated flows
  GrB_Vector e, Re;
  GrB_UnaryOp GrB_MakeFlow;
  GrB_BinaryOp GrB_InitForwardFlows, GrB_InitBackwardFlows;

  //create height vector
  GrB_Vector d;
 

  //active_set and n_active
  GrB_Vector active_set;
  GrB_Vector mask_vector;
  GrB_Index n_active;

  //semiring and vectors for y<e, struct> = R x d
  GrB_Vector y, y_dup;
  GrB_IndexUnaryOp GrB_Prune;
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
  GrB_UnaryOp GrB_extractJ, GrB_extractYJ;
  
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
 
  //relabel
  GrB_Vector d_dup;

  //update height
  GrB_BinaryOp GrB_UpdateHeight;

  //update R structure
  GrB_Matrix R_dup;
  GrB_BinaryOp GrB_UpdateFlows;

  //update e
  GrB_Vector e_dup;

  //scalars
  GrB_Scalar zero_int32;
  GrB_Scalar zero_fp32;

  //invariant
  GrB_Vector invariant;
  GrB_BinaryOp GrB_InvariantCheck;
  GrB_Scalar check;
  bool check_raw;

  //do input checks
  if(*f){
    (*f) = 0;
  }

  LG_TRY(LAGraph_CheckGraph(G, msg));
  LG_ASSERT_MSG(G->kind == LAGraph_ADJACENCY_DIRECTED, GrB_INVALID_VALUE, "LAGraph_MaxFlow requires a directed graph");

  GrB_Index ncols, nrows;
  GRB_TRY(GrB_Matrix_ncols(&ncols, G->A));
  GRB_TRY(GrB_Matrix_nrows(&nrows, G->A));
  LG_ASSERT_MSG(nrows == ncols, GrB_INVALID_VALUE, "Matrix must be square"); 
  LG_ASSERT_MSG(S < ncols && S >= 0 && T < ncols && T >= 0, GrB_INVALID_VALUE, "S and T must be a value between [0, n)");
  LG_TRY(LAGraph_Cached_EMin(G, msg));
  LG_ASSERT_MSG(G->emin > 0, GrB_INVALID_VALUE, "the edge weights (capacities) must be greater than 0");
  

  //create types for computation
  GRB_TRY(GxB_Type_new(&GrB_FlowEdge, sizeof(MF_flowEdge), "MF_flowEdge", GRB_FLOWEDGE_STR));
  GRB_TRY(GxB_Type_new(&GrB_ResultTuple, sizeof(MF_resultTuple), "MF_resultTuple", GRB_RESULTTUPLE_STR));
  GRB_TRY(GxB_Type_new(&GrB_CompareTuple, sizeof(MF_compareTuple), "MF_compareTuple", GRB_COMPARETUPLE_STR));

  //invariant check
  GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
  GRB_TRY(GxB_BinaryOp_new(&GrB_InvariantCheck, F_BINARY(MF_CheckInvariant), GrB_BOOL, GrB_INT32, GrB_ResultTuple, "MF_CheckInvariant", GRB_INV_STR));
  GRB_TRY(GrB_Scalar_new(&check, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement(check, false));
  
  //create scalars
  GRB_TRY(GrB_Scalar_new(&zero_int32, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement(zero_int32, 0));
  GRB_TRY(GrB_Scalar_new(&zero_fp32, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement(zero_fp32, 0));
  
  //create R
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidualForward, F_UNARY(MF_CreateResidualForward), GrB_FlowEdge , GrB_FP64, "MF_CreateResidualForward", GRB_CRF_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidualBackward, F_UNARY(MF_CreateResidualBackward), GrB_FlowEdge , GrB_FP64, "MF_CreateResidualBackward", GRB_CRB_STR));
  GRB_TRY(GrB_Matrix_new(&R_temp1, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_Matrix_new(&R_temp2, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_Matrix_new(&R, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_apply(R_temp1, NULL, NULL, GrB_CreateResidualForward, A, NULL));
  GRB_TRY(GrB_apply(R_temp2, NULL, NULL, GrB_CreateResidualBackward, A, GrB_DESC_T0));
  GRB_TRY(GrB_assign(R, NULL, NULL, R_temp1, GrB_ALL, n, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(R, A, NULL, R_temp2, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));

  //init R with initial saturated flows
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitForwardFlows, F_BINARY(MF_initForwardFlows), GrB_FlowEdge, GrB_FlowEdge, GrB_FlowEdge, "MF_initForwardFlows", GRB_INITFLOWF_STR));
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitBackwardFlows, F_BINARY(MF_initBackwardFlows), GrB_FlowEdge, GrB_FlowEdge, GrB_FlowEdge, "MF_initBackwardFlows", GRB_INITFLOWB_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_MakeFlow, F_UNARY(MF_MakeFlow), GrB_FlowEdge, GrB_FP64, "MF_MakeFlow", GRB_MAKEF_STR));
  GRB_TRY(GrB_Vector_new(&Re, GrB_FlowEdge, n));
  GRB_TRY(GrB_Vector_new(&e, GrB_FP64, n));
  GRB_TRY(GrB_extract(e, NULL, NULL, A, GrB_ALL, n, S, GrB_DESC_T0));
  GRB_TRY(GrB_apply(Re, NULL, NULL, GrB_MakeFlow, e, NULL));
  GRB_TRY(GxB_subassign(R, NULL, GrB_InitForwardFlows, Re, S, GrB_ALL, n, NULL));
  GRB_TRY(GxB_subassign(R, NULL, GrB_InitBackwardFlows, Re, GrB_ALL, n, S, NULL));
  
  //create and init d vector
  GRB_TRY(GrB_Vector_new(&d, GrB_INT32, n));
  GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(d, NULL, NULL, n, &S, 1, NULL));

  //extract n_active from e masking T and S then assign to e
  GRB_TRY(GrB_Vector_new(&active_set, GrB_FP64, n));
  GRB_TRY(GrB_Vector_new(&mask_vector, GrB_BOOL, n)); //keep as bool?
  GRB_TRY(GrB_assign(mask_vector, NULL, NULL, true, &T, 1, NULL));
  GRB_TRY(GrB_assign(mask_vector, NULL, NULL, true, &S, 1, NULL));
  double f_T = 0;
  GrB_Info info = GrB_Vector_extractElement(&f_T, e, T); //if value at T
  if(info == GrB_SUCCESS){
    (*f) += f_T;
  }
  GRB_TRY(GrB_select(active_set, mask_vector, NULL, GrB_VALUEGE_FP64, e, 0, GrB_DESC_RSC));
  GRB_TRY(GrB_assign(e, NULL, NULL, active_set, GrB_ALL, n, GrB_DESC_R));
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
  GRB_TRY(GxB_IndexUnaryOp_new(&GrB_Prune, (GxB_index_unary_function) MF_Prune, GrB_BOOL, GrB_ResultTuple, GrB_INT32, "MF_Prune", GRB_PRUNE_STR));

  //create utility vectors, Matrix, and ops for mapping
  GRB_TRY(GrB_Vector_new(&Jvec, GrB_INT32, n));
  GRB_TRY(GrB_Matrix_new(&map, GrB_CompareTuple, n,n));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractJ, F_UNARY(MF_extractJ), GrB_INT32, GrB_CompareTuple, "MF_extractJ", GRB_EXTRACTJ_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractYJ, F_UNARY(MF_extractYJ), GrB_INT32, GrB_ResultTuple, "MF_extractYJ", GRB_EXTRACTYJ_STR));

  //create map x e semiring
  GRB_TRY(GzB_IndexBinaryOp_new(&GrB_MxeIndexMult, F_INDEX_BINARY(MF_MxeMult), GrB_ResultTuple, GrB_CompareTuple, GrB_FP64, GrB_INT32, "MF_MxeMult", GRB_MXEMULT_STR));
  GRB_TRY(GzB_BinaryOp_new_IndexOp(&GrB_MxeMult, GrB_MxeIndexMult, theta));
  GRB_TRY(GxB_BinaryOp_new(&GrB_MxeAdd, F_BINARY(MF_MxeAdd), GrB_ResultTuple, GrB_ResultTuple, GrB_ResultTuple, "MF_MxeAdd", GRB_MXEADD_STR));
  GRB_TRY(GrB_Monoid_new_UDT(&GrB_MxeAddMonoid, GrB_MxeAdd, &id));
  GRB_TRY(GrB_Semiring_new(&GrB_MxeSemiring, GrB_MxeAddMonoid, GrB_MxeMult));

  //create flow vec
  GRB_TRY(GrB_Vector_new(&residual_vec, GrB_FP64, n));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractFlows, F_UNARY(MF_extractFlow), GrB_FP64, GrB_ResultTuple, "MF_extractFlow", GRB_EXTRACTFLOW_STR));

  GRB_TRY(GrB_Matrix_new(&delta_mat, GrB_FP64, n, n));
  GRB_TRY(GrB_Matrix_new(&delta, GrB_FP64, n, n));
  GRB_TRY(GrB_Vector_new(&delta_vec, GrB_FP64, n));

  //relable structures
  GRB_TRY(GrB_Vector_new(&d_dup, GrB_INT32, n));

  //update height binary op
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateHeight, F_BINARY(MF_updateHeight), GrB_INT32, GrB_INT32, GrB_ResultTuple, "MF_updateHeight", GRB_UPDATEHEIGHT_STR));

  //update R structure
  GRB_TRY(GrB_Matrix_new(&R_dup, GrB_FlowEdge, n, n));
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateFlows, F_BINARY(MF_updateFlow), GrB_FlowEdge, GrB_FlowEdge, GrB_FP64, "MF_updateFlow", GRB_UPDATEFLOWS_STR));

  //update e structures
  GRB_TRY(GrB_Vector_new(&e_dup, GrB_FP64, n));

  int iter = 0;
  
  while(n_active > 0){

    //BUG
    if(iter % 15 == 0 && iter > 0){
      GLOBAL_RELABEL;
    }

    //Create C arrays
    GrB_Index Jmap[LEN], Imap[LEN];
    GrB_Index Jvec_value[LEN], deltaJi[LEN];
    MF_compareTuple yd_value[LEN];
    GrB_Index Idelta[LEN], Jdelta[LEN];
    double delta_raw[LEN];


    printf("******iter: %d\n\n", iter);
    //GxB_print(e, 5);
    //GxB_print(d, 5);

    //printf("---R matrix-----\n");
    //print_flowMtx(R);
    
    //y<e, struct> = R x d
    GRB_TRY(GrB_mxv(y, e, NULL, GrB_RxdSemiring, R, d, GrB_DESC_RS));
    //printf("---y---\n\n");
    //print_resultVec(y);
    GRB_TRY(GrB_Vector_dup(&y_dup, y));
    GRB_TRY(GrB_select(y, NULL, NULL, GrB_Prune, y_dup, -1, GrB_DESC_R));

    //create yd vector of type compare tuple
    GRB_TRY(GrB_eWiseMult(yd, NULL, NULL, GrB_CreateCompareVec, y,  d, GrB_DESC_R));
    //printf("\n------yd------\n");
    //print_compareVec(yd);

    //create map matrix from yd
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractJ, yd, GrB_DESC_R));
    //GxB_print(Jvec, 5);
    GrB_Index JVec_n, yd_n;
    GRB_TRY(GrB_Vector_nvals(&JVec_n, Jvec));
    GRB_TRY(GrB_Vector_nvals(&yd_n, yd));
    GRB_TRY(GrB_Vector_extractTuples(Jmap, Jvec_value, &JVec_n, Jvec));
    GRB_TRY(GrB_Vector_extractTuples(Imap, (void*)yd_value, &yd_n, yd));
    GRB_TRY(GrB_Matrix_build(map, Imap, Jvec_value, (void*)yd_value, yd_n, GxB_IGNORE_DUP));
    //GxB_print(map, 5);
    
    //make e dense for map computation
    GRB_TRY(GrB_assign(e, e, NULL, 0, GrB_ALL, n, GrB_DESC_SC));
    //printf("***begin map computation***\n\n");
    //GxB_print(e, 5);

    //y = map x e
    GRB_TRY(GrB_mxv(y_dup, NULL, NULL, GrB_MxeSemiring, map, e, GrB_DESC_R));
    //printf("******MAP***********\n\n");
    //print_MapMtx(map);
    //printf("\n");
    //printf("----y-prePrune----\n\n");
    //print_resultVec(y_dup);
    GRB_TRY(GrB_select(y, NULL, NULL, GrB_Prune, y_dup, -1, GrB_DESC_R));
    //printf("----y-postPrune----\n\n");
    //print_resultVec(y);

    //relable, update heights
    // add alpha and beta scalars
    GRB_TRY(GrB_Vector_dup(&d_dup, d));
    GRB_TRY(GrB_eWiseMult(d, y, NULL, GrB_UpdateHeight, d_dup, y, GrB_DESC_S));
    //GxB_print(d, 5);
    //assert correct labels
    GRB_TRY(GrB_eWiseMult(invariant, y, NULL, GrB_InvariantCheck, d, y, GrB_DESC_RS));
    GRB_TRY(GrB_reduce(check, NULL, GrB_LAND_MONOID_BOOL, invariant, GrB_DESC_R));
    GRB_TRY(GrB_Scalar_extractElement(&check_raw, check));
    //GxB_print(d, 5);
    //printf("\n");
    //print_resultVec(y);
    //GxB_print(invariant, 5);
    LG_ASSERT_MSG(check_raw == true, GrB_PANIC, "The invariant is not upheld, the algorithm is wrong!!");
    
    //GxB_print(d, 5);

    //extract residual flows from y
    GRB_TRY(GrB_apply(residual_vec, NULL, NULL, GrB_extractFlows, y, GrB_DESC_R));
    //GxB_print(residual_vec, 5);

    //.min(flow_vec and e)
    GRB_TRY(GrB_eWiseMult(delta_vec, NULL, NULL, GrB_MIN_FP64, residual_vec, e, GrB_DESC_R));
    //GxB_print(delta_vec, 5);

    //create delta matrix from delta vector
    GrB_Index delta_vec_n, delta_jn;
    GRB_TRY(GrB_Vector_nvals(&delta_vec_n, delta_vec));
    GRB_TRY(GrB_Vector_extractTuples(Idelta, delta_raw, &delta_vec_n, delta_vec));
    GRB_TRY(GrB_Vector_nvals(&delta_jn, y));
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractYJ, y, GrB_DESC_R));
    GRB_TRY(GrB_Vector_extractTuples(deltaJi, Jdelta, &delta_jn, Jvec));
    GRB_TRY(GrB_Matrix_build(delta, Idelta, Jdelta, delta_raw, delta_vec_n, GxB_IGNORE_DUP));
    //GxB_print(delta, 5);

    //make delta anti-symmetric
    GRB_TRY(GxB_eWiseUnion(delta_mat, NULL, NULL, GrB_MINUS_FP64, delta, zero_fp32, delta, zero_fp32, GrB_DESC_RT1));
    //GxB_print(delta_mat, 5);

    //update R
    GRB_TRY(GrB_Matrix_dup(&R_dup, R));
    GRB_TRY(GrB_eWiseMult(R, delta_mat, NULL, GrB_UpdateFlows, R_dup, delta_mat, GrB_DESC_S));

    //printf("---New R---\n");
    //print_flowMtx(R);

    //reduce delat_mat to delta_vec
    GRB_TRY(GrB_reduce(delta_vec, NULL, NULL, GrB_PLUS_FP64, delta_mat, GrB_DESC_RT0));
    //GxB_print(delta_vec, 5);

    //add to e
    //add alpha and beta scalars
    GRB_TRY(GrB_Vector_dup(&e_dup, e));
    GRB_TRY(GxB_eWiseUnion(e, NULL, NULL, GrB_PLUS_FP64, e_dup, zero_fp32, delta_vec, zero_fp32, NULL));
    //GxB_print(e, 5);

    //TO-DO:
    //backwards BFS periodically, if can't find src: give up
    //add checks for d vector 
    //
    
    //extract active nodes
    GRB_TRY(GrB_Vector_extractElement(&f_T, e, T));
    (*f) += f_T;
    GRB_TRY(GrB_select(active_set, mask_vector, NULL, GrB_VALUEGT_FP64, e, 0, GrB_DESC_RSC));
    GRB_TRY(GrB_assign(e, NULL, NULL, active_set, GrB_ALL, n, GrB_DESC_R));
    GRB_TRY(GrB_Vector_nvals(&n_active, active_set));
    //GxB_print(active_set, 5);
    //printf("max flow in alg iter is: %f\n", *f);

    //clear map and delta
    GRB_TRY(GrB_Matrix_clear(map));
    GRB_TRY(GrB_Matrix_clear(delta));

    ++iter;
    
  }

  //print_flowMtx(R);
  //GxB_print(d, 5);
  
  LG_FREE_ALL;
  return GrB_SUCCESS;
}
