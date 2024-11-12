
#define LG_FREE_WORK                        \
{                                           \
    /* free any workspace used here */      \
    GrB_free (&W) ;                         \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (&Y) ;                         \
    /* take any other corrective action */  \
}

#include <LG_internal.h>
#include <LAGraphX.h>
#include <LAGraph.h>

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

void Rxd_MultBOp(MF_result_tuple *y, GrB_Flow_Edge *R, int *d, int i, int j) {
  float r = R->capacity - R->flow;
  if(r > 0){
    y->residual = r;
    y->d = *d;
    y->j = j;
  }
  else{
    y->residual = 0;
    y->d = INT32_MAX;
    y->j = -1;
  }
}

void Rxd_AddMonoid(MF_result_tuple * z, MF_result_tuple * x, MF_result_tuple * y) {
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
int LAGraph_MaxFlow(LAGraph_Graph *G, GrB_Index S, GrB_Index T, GrB_Matrix *R, int * f, char *msg){
  
  
  return 0;
}
