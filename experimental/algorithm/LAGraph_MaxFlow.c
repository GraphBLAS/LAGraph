//************************HOUSE RULES *******************
// user defined typedefs are MF for Max Flow
// when registered with GraphBLAS, replace MF with GrB
// macros are all caps
// op params will be z, y, x accordingly

#include <LAGraphX.h>
#include "LG_internal.h"
#include <LAGraph.h>

// FIXME: replace int with int64, and INT32 with INT64.
// Consider two sets of data types and ops: 32/64 bit integers ...

//------------------------------------------------------------------------------
// LG_augment_maxflow
//------------------------------------------------------------------------------

// FIXME: describe me

#undef  LG_FREE_ALL
#define LG_FREE_ALL ;

static GrB_Info LG_augment_maxflow
(
    double *f,              // maxflow
    GrB_Vector e,
    GrB_Index T,            // sink node
    GrB_Vector mask_vector,
    GrB_Vector active_set,
    GrB_Index *n_active,
    GrB_Index n,
    char *msg
)
{
    // f_T = e (T)
    double f_T = 0;
    GrB_Info info = GrB_Vector_extractElement(&f_T, e, T); //if value at T
    GRB_TRY (info) ;
    if (info == GrB_SUCCESS)
    {
        // e(T) is present
        (*f) += f_T;
    }
    GRB_TRY(GrB_select(active_set, mask_vector, NULL, GrB_VALUEGT_FP64, e, 0, GrB_DESC_RSC));
    // e = active_set : FIXME do we need this?
    GRB_TRY(GrB_assign(e, NULL, NULL, active_set, GrB_ALL, n, GrB_DESC_R));
    GRB_TRY(GrB_Vector_nvals(n_active, active_set));
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
    GrB_free(&delta);                                                          \
    GrB_free(&delta_vec);                                                      \
    GrB_free(&delta_mat);                                                      \
    GrB_free(&active_set);                                                     \
    GrB_free(&map);                                                            \
    GrB_free(&y);                                                              \
    GrB_free(&yd);                                                             \
    GrB_free(&mask_vector);                                                    \
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
    GrB_free(&zero_fp64);                                                      \
    GrB_free(&Re);                                                             \
    GrB_free(&invariant);                                                      \
    GrB_free(&GrB_InvariantCheck);                                             \
    GrB_free(&check);                                                          \
    GrB_free(&GrB_extractYJ);                                                  \
    GrB_free(&extract_desc);                                                   \
    GrB_free(&residual_vec);						\
    GrB_free(&GrB_MakeFlow);						\
    GrB_free(&GrB_GetResidual);				\
  }


#define LG_FREE_ALL \
{ \
  LG_FREE_WORK; \
}

//#define LEN INT16_MAX

//casting for unary ops
#define F_UNARY(f) ((void (*)(void *, const void *))f)

// casting for index binary ops
#define F_INDEX_BINARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void *, GrB_Index, GrB_Index, const void *)) f)

// casting for index unary op
//#define F_INDEX_UNARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void*)) f)

// casting for binary op
#define F_BINARY(f) ((void (*)(void *, const void *, const void *)) f)

//custom types
typedef struct{
  double flow;
  double capacity;
} MF_flowEdge;

#define GRB_FLOWEDGE_STR "typedef struct{"\
  "double flow;"\
  "double capacity;"\
"} MF_flowEdge;"


typedef struct{
  double residual;
  GrB_Index j;
  int64_t d;
} MF_resultTuple;

#define GRB_RESULTTUPLE_STR "typedef struct{"\
  "double residual;"\
  "GrB_Index j;"\
  "int64_t d;"\
"} MF_resultTuple;"

typedef struct{
  double residual;
  int64_t di;
  int64_t y_dmin;
  GrB_Index j;
} MF_compareTuple;

#define GRB_COMPARETUPLE_STR                                                   \
  "typedef struct{"                                                            \
  "double residual;"                                                            \
  "int64_t di;"                                                                    \
  "int64_t y_dmin;"                                                                \
  "GrB_Index j;"                                                               \
  "} MF_compareTuple;"

void MF_CreateResidualForward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = (*y);
}

#define GRB_CRF_STR "void MF_CreateResidualForward(MF_flowEdge *z, const double *y) {"\
  "z->flow = 0;"\
  "z->capacity = (*y);"\
"}"

void MF_CreateResidualBackward(MF_flowEdge *z, const double *y) {
  z->flow = 0;
  z->capacity = 0;
}

#define GRB_CRB_STR "void MF_CreateResidualBackward(MF_flowEdge *z, const double *y)"\
 "{"\
  "z->flow " \
  "= 0;\nz->capacity = 0;"\
"}"


void MF_RxdMult(MF_resultTuple *z, const MF_flowEdge *y, GrB_Index iy, GrB_Index jy, const int64_t *x, GrB_Index ix, GrB_Index jx, const int64_t* theta) {
  double r = y->capacity - y->flow;
  if(r > 0){
    z->residual = r;
    z->d = *x;
    z->j = jy;
  }
  else if(r==0){
    z->residual = 0;
    z->d = INT64_MAX;
    z->j = -1;
  }
}

#define GRB_RXDMULT_STR "void MF_RxdMult(MF_resultTuple *z, const MF_flowEdge *y, GrB_Index iy, GrB_Index jy, const int64_t *x, GrB_Index ix, GrB_Index jx, const int64_t* theta) {"\
  "double r = y->capacity - y->flow;"\
  "if(r > 0){"\
    "z->residual = r;"\
    "z->d = *x;"\
    "z->j = jy;"\
  "}"\
  "else{"\
    "z->residual = 0;"\
    "z->d = INT64_MAX;"\
    "z->j = -1;"\
  "}"\
"}"


void MF_RxdAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x) {
  if(y->d < x->d){
    (*z) = (*y) ;
  }
  else if(y->d > x->d){
    (*z) = (*x) ;
  }
  else{
    if(y->residual > x->residual){
      (*z) = (*y) ;
    }
    else if(y->residual < x->residual){
      (*z) = (*x) ;
    }
    else{
      if(y->j > x->j){
	(*z) = (*y);
      }
      else{
	(*z) = (*x) ;
      }
    }
  }
}

#define GRB_RXDADD_STR "void MF_RxdAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x) {"\
  "if(y->d < x->d){"\
    "(*z) = (*y) ;"\
  "}"\
  "else if(y->d > x->d){"\
    "(*z) = (*x) ;"\
  "}"\
  "else{"\
    "if(y->residual > x->residual){"\
      "(*z) = (*y) ;"\
    "}"\
    "else if(y->residual < x->residual){"\
      "(*z) = (*x) ;"\
    "}"\
    "else{"\
      "if(y->j > x->j){"\
	"(*z) = (*y) ;"\
      "}"\
      "else{"\
	"(*z) = (*x) ;"\
      "}"\
    "}"\
  "}"\
"}"


void MF_extractFlow(double *z, const MF_resultTuple *y) { (*z) = y->residual; }

#define GRB_EXTRACTFLOW_STR "void MF_extractFlow(double *z, const MF_resultTuple *y) {\n(*z) = y->residual;\n}"

void MF_updateFlow(MF_flowEdge *z, const MF_flowEdge *y, const double *x) {
  z->capacity = y->capacity;
  z->flow = y->flow + (*x);
}

#define GRB_UPDATEFLOWS_STR "void MF_updateFlow(MF_flowEdge *z, const MF_flowEdge *y, const double *x) {"\
  "z->capacity = y->capacity;"\
  "z->flow = y->flow + (*x);"\
"}"


void MF_updateHeight(int64_t *z, const int64_t *y, const MF_resultTuple *x) {
  if((*y) < x->d+1){
    (*z) = x->d + 1;
  }
  else if ((*y) == x->d+1){
    (*z) = (*y);
  }
}

#define GRB_UPDATEHEIGHT_STR "void MF_updateHeight(int64_t *z, const int64_t *y, const MF_resultTuple *x) {"\
  "if((*y) < x->d+1){"\
    "(*z) = x->d + 1;"\
  "}"\
  "else if ((*y) == x->d+1){"\
    "(*z) = (*y);"\
  "}"\
"}"


void MF_extractJ(int64_t *z, const MF_compareTuple *y) { (*z) = y->j; }

#define GRB_EXTRACTJ_STR "void MF_extractJ(int64_t *z, const MF_compareTuple *y) {(*z) = y->j;}"

void MF_extractYJ(int64_t *z, const MF_resultTuple *y) {
  (*z) = y->j;
}

#define GRB_EXTRACTYJ_STR "void MF_extractYJ(int64_t *z, const MF_resultTuple *y) {" \
  "(*z) = y->j;" \
"}"


void MF_initForwardFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){
  z->flow = x->flow + y->flow;
  z->capacity = y->capacity;
}

#define GRB_INITFLOWF_STR "void MF_initForwardFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){"\
  "z->flow = x->flow + y->flow;"\
  "z->capacity = y->capacity;"\
"}"


void MF_initBackwardFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){
  z->flow = y->flow - x->flow;
  z->capacity = y->capacity;
}

#define GRB_INITFLOWB_STR "void MF_initBackwardFlows(MF_flowEdge * z, const MF_flowEdge * y, const MF_flowEdge * x){"\
  "z->flow = y->flow - x->flow;"\
  "z->capacity = y->capacity;"\
"}"


void MF_MxeMult(MF_resultTuple * z, const MF_compareTuple * y, GrB_Index iy, GrB_Index jy, const double * x, GrB_Index ix, GrB_Index jx, const int64_t* theta){
  if(y->di == y->y_dmin && (*x) > 0){ //check this
    if(iy < jy){
      z->d = y->y_dmin;
      z->residual = y->residual;
      z->j = y->j;
    } //add else to populate with empty tuple, prune after.
    else{
      z->d = INT64_MAX;
      z->j = -1;
      z->residual = 0;
    }
  }
  else if(y->di == y->y_dmin - 1 && (*x) > 0){
    z->d = INT64_MAX;
    z->j = -1;
    z->residual = 0;
  }
  else if(y->di < y->y_dmin-1 || y->di == y->y_dmin+1 || y->di == y->y_dmin-1 || y->di == y->y_dmin){
    z->d = y->y_dmin;
    z->residual = y->residual;
    z->j = y->j;
  }
  else{ //change later to signify the removal of the node from the active set since flow cannot be pushed anywhere.
    z->d = INT64_MAX;
    z->residual = 0;
    z->j = -1;
  }
}

#define GRB_MXEMULT_STR "void MF_MxeMult(MF_resultTuple * z, const MF_compareTuple * y, GrB_Index iy, GrB_Index jy, const double * x, GrB_Index ix, GrB_Index jx, const int64_t* theta){" \
  "if(y->di == y->y_dmin && (*x) > 0){" \
    "if(iy < jy){" \
      "z->d = y->y_dmin;" \
      "z->residual = y->residual;" \
      "z->j = y->j;" \
    "}" \
    "else{"  \
      "z->d = INT64_MAX;" \
      "z->j = -1;" \
      "z->residual = 0;" \
    "}" \
  "}" \
  "else if(y->di == y->y_dmin - 1 && (*x) > 0){" \
    "z->d = INT64_MAX;" \
    "z->j = -1;" \
    "z->residual = 0;" \
  "}" \
  "else if(y->di < y->y_dmin-1 || y->di == y->y_dmin+1 || y->di == y->y_dmin-1 || y->di == y->y_dmin){" \
    "z->d = y->y_dmin;" \
    "z->residual = y->residual;" \
    "z->j = y->j;" \
  "}" \
  "else{" \
    "z->d = INT64_MAX;" \
    "z->residual = 0;" \
    "z->j = -1;" \
  "}" \
"}"


void MF_MxeAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x){
  if(x != NULL){
    (*z) = (*x) ;
  }
  else{
    (*z) = (*y) ;
  }
}

#define GRB_MXEADD_STR "void MF_MxeAdd(MF_resultTuple * z, const MF_resultTuple * y, const MF_resultTuple * x){"\
  "if(x != NULL){"\
    "(*z) = (*x) ;"\
  "}"\
  "else{"\
    "(*z) = (*y) ;"\
  "}"\
"}"


void MF_CreateCompareVec(MF_compareTuple *z, const MF_resultTuple *y, const int64_t *x) {
  z->di = (*x);
  z->j = y->j;
  z->residual = y->residual;
  z->y_dmin = y->d;
}

#define GRB_CREATECOMPVEC_STR "void MF_CreateCompareVec(MF_compareTuple *z, const MF_resultTuple *y, const int64_t *x) {"\
  "z->di = (*x);"\
  "z->j = y->j;"\
  "z->residual = y->residual;"\
  "z->y_dmin = y->d;"\
"}"


void MF_Prune(bool * z, const MF_resultTuple * y, GrB_Index iy, GrB_Index jy, const int64_t * theta){
  if(y->j != *theta){
    *z = true;
  }
  else{
    *z = false;
  }
}

#define GRB_PRUNE_STR "void MF_Prune(bool * z, const MF_resultTuple * y, GrB_Index iy, GrB_Index jy, const int64_t * theta){"\
  "if(y->j != *theta){"\
    "*z = true;" \
  "}" \
  "else{" \
    "*z = false;" \
  "}" \
"}"

void MF_MakeFlow(MF_flowEdge * z, const double * y){
  z->capacity = 0;
  z->flow = (*y);
}

#define GRB_MAKEF_STR "void MF_MakeFlow(MF_flowEdge * z, const double * y){"\
  "z->capacity = 0;"\
  "z->flow = (*y);"\
"}"

// FIXME: fix GraphBLAS so it can print user-defined types
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
    printf("(%ld, %ld)         (height: %ld, y.dmin: %ld, J: %ld, residual: %lf) \n", i, j, e.di, e.y_dmin, e.j, e.residual);
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
    printf("(%ld, 0)         (height: %ld, J: %ld, residual: %f) \n", i, e.d, e.j, e.residual);
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
    printf("(%ld, 0)         (height: %ld, y.dmin: %ld, J: %ld, residual: %lf) \n", i, e.di, e.y_dmin, e.j, e.residual);
    info = GxB_Vector_Iterator_next(iter);
  }
  GrB_free(&iter);
}


void MF_CheckInvariant(bool *z, const int *y, const MF_resultTuple *x) {
  (*z) = ((*y) == x->d+1);
}

#define GRB_INV_STR "void MF_CheckInvariant(bool *z, const int *y, const MF_resultTuple *x) {"\
  "(*z) = ((*y) == x->d+1);"\
"}"


void MF_getResidual(double * z, const MF_flowEdge * y){
  *z = y->capacity - y->flow;
}

#define GRB_GETRES_STR "void MF_getResidual(double * z, const MF_flowEdge * y){" \
  "*z = y->capacity - y->flow;" \
"}"
  
//------------------------------------------------------------------------------
// LAGraph_MaxFlow
//------------------------------------------------------------------------------

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
  GrB_Type GrB_FlowEdge = NULL ;
  GrB_Type GrB_ResultTuple = NULL ;
  GrB_Type GrB_CompareTuple = NULL ;

  GrB_Vector lvl = NULL ;
  GrB_UnaryOp GrB_GetResidual = NULL ;
  GrB_Matrix res_mat = NULL, res_matT = NULL ;
  LAGraph_Graph res_graph = NULL ;
    

  //to create R
  GrB_UnaryOp GrB_CreateResidualForward = NULL , GrB_CreateResidualBackward = NULL ;
  GrB_Matrix A = G->A;  /* FIXME, move below */
  GrB_Index n;
  GrB_Matrix R = NULL ;
  GrB_Matrix_nrows(&n, A);

  //to init R with initial saturated flows
  GrB_Vector e = NULL, Re = NULL ;
  GrB_UnaryOp GrB_MakeFlow = NULL ;
  GrB_BinaryOp GrB_InitForwardFlows = NULL, GrB_InitBackwardFlows = NULL ;

  //create height vector
  GrB_Vector d = NULL ;
 

  //active_set and n_active
  GrB_Vector active_set = NULL ;
  GrB_Vector mask_vector = NULL ;
  GrB_Index n_active ;

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
 
  //delta structures
  GrB_Vector delta_vec = NULL ;
  GrB_Matrix delta = NULL , delta_mat = NULL ;
 
  //update height
  GrB_BinaryOp GrB_UpdateHeight = NULL ;

  //update R structure
  GrB_BinaryOp GrB_UpdateFlows = NULL ;

  //scalars
  GrB_Scalar zero_fp64 = NULL ;

  //invariant
  GrB_Vector invariant = NULL ;
  GrB_BinaryOp GrB_InvariantCheck = NULL ;
  GrB_Scalar check = NULL ;
  bool check_raw;

  //descriptor and matrix building
  GrB_Descriptor extract_desc = NULL ;

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
  LG_ASSERT_MSG(G->emin > 0, GrB_INVALID_VALUE, "the edge weights (capacities) must be greater than 0");
  
  //create types for computation
  GRB_TRY(GxB_Type_new(&GrB_FlowEdge, sizeof(MF_flowEdge), "MF_flowEdge", GRB_FLOWEDGE_STR));
  GRB_TRY(GxB_Type_new(&GrB_ResultTuple, sizeof(MF_resultTuple), "MF_resultTuple", GRB_RESULTTUPLE_STR));
  GRB_TRY(GxB_Type_new(&GrB_CompareTuple, sizeof(MF_compareTuple), "MF_compareTuple", GRB_COMPARETUPLE_STR));

  //global relabel operations
  GRB_TRY(GxB_UnaryOp_new(&GrB_GetResidual, F_UNARY(MF_getResidual), GrB_FP64, GrB_FlowEdge, "MF_getResidual", GRB_GETRES_STR));
  
  
  //invariant check
  GRB_TRY(GrB_Vector_new(&invariant, GrB_BOOL, n));
  GRB_TRY(GxB_BinaryOp_new(&GrB_InvariantCheck, F_BINARY(MF_CheckInvariant), GrB_BOOL, GrB_INT64, GrB_ResultTuple, "MF_CheckInvariant", GRB_INV_STR));
  GRB_TRY(GrB_Scalar_new(&check, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement(check, false));
  
  //create scalars
  GRB_TRY(GrB_Scalar_new(&zero_fp64, GrB_FP64));
  GRB_TRY(GrB_Scalar_setElement(zero_fp64, 0));
  
  //create R
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidualForward, F_UNARY(MF_CreateResidualForward), GrB_FlowEdge , GrB_FP64, "MF_CreateResidualForward", GRB_CRF_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_CreateResidualBackward, F_UNARY(MF_CreateResidualBackward), GrB_FlowEdge , GrB_FP64, "MF_CreateResidualBackward", GRB_CRB_STR));
  GRB_TRY(GrB_Matrix_new(&R, GrB_FlowEdge, n, n));
  GRB_TRY(GrB_apply(R, NULL, NULL, GrB_CreateResidualForward, A, NULL));
  //FIXME: rename to LAGr_MaxFlow, and use G->AT here:
  GRB_TRY(GrB_apply(R, A, NULL, GrB_CreateResidualBackward, G->AT, GrB_DESC_SC));

  //init R with initial saturated flows
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitForwardFlows, F_BINARY(MF_initForwardFlows), GrB_FlowEdge, GrB_FlowEdge, GrB_FlowEdge, "MF_initForwardFlows", GRB_INITFLOWF_STR));
  GRB_TRY(GxB_BinaryOp_new(&GrB_InitBackwardFlows, F_BINARY(MF_initBackwardFlows), GrB_FlowEdge, GrB_FlowEdge, GrB_FlowEdge, "MF_initBackwardFlows", GRB_INITFLOWB_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_MakeFlow, F_UNARY(MF_MakeFlow), GrB_FlowEdge, GrB_FP64, "MF_MakeFlow", GRB_MAKEF_STR));
  GRB_TRY(GrB_Vector_new(&Re, GrB_FlowEdge, n));
  GRB_TRY(GrB_Vector_new(&e, GrB_FP64, n));
  GRB_TRY(GrB_extract(e, NULL, NULL, A, GrB_ALL, n, S, GrB_DESC_T0));
  GRB_TRY(GrB_apply(Re, NULL, NULL, GrB_MakeFlow, e, NULL));
  GRB_TRY(GrB_assign(R, NULL, GrB_InitForwardFlows, Re, S, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(R, NULL, GrB_InitBackwardFlows, Re, GrB_ALL, n, S, NULL));
  
  //create and init d vector
  GRB_TRY(GrB_Vector_new(&d, GrB_INT64, n));
  GRB_TRY(GrB_assign(d, NULL, NULL, 0, GrB_ALL, n, NULL));
  GRB_TRY(GrB_assign(d, NULL, NULL, n, &S, 1, NULL));

  //extract n_active from e masking T and S then assign to e
  GRB_TRY(GrB_Vector_new(&active_set, GrB_FP64, n));
  GRB_TRY(GrB_Vector_new(&mask_vector, GrB_BOOL, n)); //keep as bool?
  GRB_TRY (GrB_Vector_setElement (mask_vector, true, T)) ;
  GRB_TRY (GrB_Vector_setElement (mask_vector, true, S)) ;

  // augment maxflow if the edge (S,T) exists
  LG_TRY (LG_augment_maxflow (f, e, T, mask_vector, active_set, &n_active, n, msg)) ;

  //create semiring and vectors for y<e, struct> = R x d
  GRB_TRY(GrB_Scalar_new(&theta, GrB_INT32));
  GRB_TRY(GrB_Scalar_setElement_INT64(theta, 0));
  GRB_TRY(GrB_Vector_new(&y, GrB_ResultTuple, n));
  GRB_TRY(GxB_IndexBinaryOp_new(&GrB_RxdIndexMult, F_INDEX_BINARY(MF_RxdMult), GrB_ResultTuple, GrB_FlowEdge, GrB_INT64, GrB_INT64, "MF_RxdMult", GRB_RXDMULT_STR));
  GRB_TRY(GxB_BinaryOp_new_IndexOp(&GrB_RxdMult, GrB_RxdIndexMult, theta));
  GRB_TRY(GxB_BinaryOp_new(&GrB_RxdAdd, F_BINARY(MF_RxdAdd), GrB_ResultTuple, GrB_ResultTuple, GrB_ResultTuple, "MF_RxdAdd", GRB_RXDADD_STR));
  MF_resultTuple id = {.d = INT64_MAX, .j = -1, .residual = 0};

  GRB_TRY(GrB_Monoid_new_UDT(&GrB_RxdAddMonoid, GrB_RxdAdd, &id));
  GRB_TRY(GrB_Semiring_new(&GrB_RxdSemiring, GrB_RxdAddMonoid, GrB_RxdMult));

  //create binary op and yd
  GRB_TRY(GrB_Vector_new(&yd, GrB_CompareTuple, n));
  GRB_TRY(GxB_BinaryOp_new(&GrB_CreateCompareVec, F_BINARY(MF_CreateCompareVec), GrB_CompareTuple, GrB_ResultTuple, GrB_INT64, "MF_CreateCompareVec", GRB_CREATECOMPVEC_STR));
  GRB_TRY(GxB_IndexUnaryOp_new(&GrB_Prune, (GxB_index_unary_function) MF_Prune, GrB_BOOL, GrB_ResultTuple, GrB_INT64, "MF_Prune", GRB_PRUNE_STR));

  //create utility vectors, Matrix, and ops for mapping
  GrB_Type JType = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32;
  GRB_TRY(GrB_Vector_new(&Jvec, JType, n));
  GRB_TRY(GrB_Matrix_new(&map, GrB_CompareTuple, n,n));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractJ, F_UNARY(MF_extractJ), GrB_INT64, GrB_CompareTuple, "MF_extractJ", GRB_EXTRACTJ_STR));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractYJ, F_UNARY(MF_extractYJ), GrB_INT64, GrB_ResultTuple, "MF_extractYJ", GRB_EXTRACTYJ_STR));

  //create map x e semiring
  GRB_TRY(GxB_IndexBinaryOp_new(&GrB_MxeIndexMult, F_INDEX_BINARY(MF_MxeMult), GrB_ResultTuple, GrB_CompareTuple, GrB_FP64, GrB_INT64, "MF_MxeMult", GRB_MXEMULT_STR));
  GRB_TRY(GxB_BinaryOp_new_IndexOp(&GrB_MxeMult, GrB_MxeIndexMult, theta));
  GRB_TRY(GxB_BinaryOp_new(&GrB_MxeAdd, F_BINARY(MF_MxeAdd), GrB_ResultTuple, GrB_ResultTuple, GrB_ResultTuple, "MF_MxeAdd", GRB_MXEADD_STR));
  GRB_TRY(GrB_Monoid_new_UDT(&GrB_MxeAddMonoid, GrB_MxeAdd, &id));
  GRB_TRY(GrB_Semiring_new(&GrB_MxeSemiring, GrB_MxeAddMonoid, GrB_MxeMult));

  //create flow vec
  GRB_TRY(GrB_Vector_new(&residual_vec, GrB_FP64, n));
  GRB_TRY(GxB_UnaryOp_new(&GrB_extractFlows, F_UNARY(MF_extractFlow), GrB_FP64, GrB_ResultTuple, "MF_extractFlow", GRB_EXTRACTFLOW_STR));

  GRB_TRY(GrB_Matrix_new(&delta_mat, GrB_FP64, n, n));
  GRB_TRY(GrB_Matrix_new(&delta, GrB_FP64, n, n));
  GRB_TRY(GrB_Vector_new(&delta_vec, GrB_FP64, n));

  //update height binary op
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateHeight, F_BINARY(MF_updateHeight), GrB_INT64, GrB_INT64, GrB_ResultTuple, "MF_updateHeight", GRB_UPDATEHEIGHT_STR));

  //update R structure
  GRB_TRY(GxB_BinaryOp_new(&GrB_UpdateFlows, F_BINARY(MF_updateFlow), GrB_FlowEdge, GrB_FlowEdge, GrB_FP64, "MF_updateFlow", GRB_UPDATEFLOWS_STR));

  int64_t iter = 0;

  //Create extract arrays
  GRB_TRY(GrB_Descriptor_new(&extract_desc));
  GRB_TRY(GrB_set(extract_desc, GxB_USE_INDICES, GxB_ROWINDEX_LIST)); 
  
  while(n_active > 0){

    // global relabeling, for the first iteration, and every 12 iterations
    // after that
    if(iter % 12 == 0)
    {
      GRB_TRY(GrB_Matrix_new(&res_matT, GrB_FP64, n, n));
      GRB_TRY(GrB_Matrix_new(&res_mat, GrB_FP64, n, n));
      GRB_TRY(GrB_apply(res_mat, NULL, NULL, GrB_GetResidual, R, NULL)) ;
      GRB_TRY(GrB_select(res_mat, NULL, NULL, GrB_VALUEGT_FP64, res_mat, 0, NULL)) ;
      GRB_TRY(GrB_transpose(res_matT, NULL, NULL, res_mat, NULL));
      LG_TRY(LAGraph_New(&res_graph, &res_matT, LAGraph_ADJACENCY_DIRECTED, msg));
      res_graph->AT = res_mat;
      res_mat = NULL ;
      LG_TRY(LAGraph_Cached_OutDegree(res_graph, msg));
      LG_TRY(LAGr_BreadthFirstSearch(&lvl, NULL, res_graph, T, msg));
      GRB_TRY(GrB_assign(d, mask_vector, NULL, lvl, GrB_ALL, n, GrB_DESC_SC));
      GRB_TRY(GrB_assign(d, lvl, NULL, n, GrB_ALL, n, GrB_DESC_SC));
      GRB_TRY(GrB_assign(e, lvl, NULL, -1, GrB_ALL, n, GrB_DESC_SC));
      GRB_TRY(GrB_select(e, NULL, NULL, GrB_VALUEGT_FP64, e, -1, NULL));
      GrB_free(&lvl);
      LG_TRY(LAGraph_Delete(&res_graph, msg));
      GRB_TRY(GrB_Vector_nvals(&n_active, e));
      if(n_active == 0){
	printf("exited early!\n");  // FIXME remove printfs
	break;
      }
    }

    printf("******iter: %ld\n\n", iter); 
    
    GRB_TRY(GrB_mxv(y, e, NULL, GrB_RxdSemiring, R, d, GrB_DESC_RS));
    GRB_TRY(GrB_select(y, NULL, NULL, GrB_Prune, y, -1, NULL));

    //create yd vector of type compare tuple
    GRB_TRY(GrB_eWiseMult(yd, NULL, NULL, GrB_CreateCompareVec, y,  d, NULL));

    //create map matrix from yd
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractJ, yd, NULL));
    GRB_TRY(GrB_Matrix_clear(map));
    GRB_TRY(GrB_Matrix_build(map, yd, Jvec, yd, GxB_IGNORE_DUP, extract_desc));
    
    //make e dense for map computation
    GRB_TRY(GrB_assign(e, e, NULL, 0, GrB_ALL, n, GrB_DESC_SC));

    //y = map x e
    GRB_TRY(GrB_mxv(y, NULL, NULL, GrB_MxeSemiring, map, e, NULL));
    GRB_TRY(GrB_select(y, NULL, NULL, GrB_Prune, y, -1, NULL));
   
    //relable, update heights
    GRB_TRY(GrB_eWiseMult(d, y, NULL, GrB_UpdateHeight, d, y, GrB_DESC_S));

    #ifdef DBG
    //assert correct labels
        GRB_TRY(GrB_eWiseMult(invariant, y, NULL, GrB_InvariantCheck, d, y, GrB_DESC_RS));
	GRB_TRY(GrB_reduce(check, NULL, GrB_LAND_MONOID_BOOL, invariant, NULL));
	GRB_TRY(GrB_Scalar_extractElement(&check_raw, check));
	ASSERT(check_raw == true);
    #endif

    //extract residual flows from y
    GRB_TRY(GrB_apply(residual_vec, NULL, NULL, GrB_extractFlows, y, NULL));

    //.min(flow_vec and e)
    GRB_TRY(GrB_eWiseMult(delta_vec, NULL, NULL, GrB_MIN_FP64, residual_vec, e, NULL));
    GRB_TRY(GrB_apply(Jvec, NULL, NULL, GrB_extractYJ, y, NULL));
    GRB_TRY(GrB_Matrix_clear(delta));
    GRB_TRY(GxB_Matrix_build_Vector(delta, delta_vec, Jvec, delta_vec, GxB_IGNORE_DUP, extract_desc));

    //make delta anti-symmetric
    GRB_TRY(GxB_eWiseUnion(delta_mat, NULL, NULL, GrB_MINUS_FP64, delta, zero_fp64, delta, zero_fp64, GrB_DESC_T1));

    //update R
    GRB_TRY(GrB_eWiseMult(R, delta_mat, NULL, GrB_UpdateFlows, R, delta_mat, GrB_DESC_S));

    //reduce delta_mat to delta_vec
    GRB_TRY(GrB_reduce(delta_vec, NULL, NULL, GrB_PLUS_FP64, delta_mat, GrB_DESC_T0));

    //add to e
    GRB_TRY(GrB_assign(e, delta_vec, GrB_PLUS_FP64, delta_vec, GrB_ALL, n, GrB_DESC_S));
    
    // augment maxflow for all active nodes
    LG_TRY (LG_augment_maxflow (f, e, T, mask_vector, active_set, &n_active, n, msg)) ;

    ++iter;
    
  }

  //print_flowMtx(R);
  printf("DBG: number of active = %ld\n", n_active);
  LG_FREE_ALL;
  return GrB_SUCCESS;
}
