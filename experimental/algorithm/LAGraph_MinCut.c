#include <LAGraph.h>
#include "LG_internal.h"
#include <LAGraph.h>

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK				\
{						\
 GrB_free(&S_diag);				\
 GrB_free(&S_bar_diag);				\
 LAGraph_Delete(&G, msg);			\
}

#define LG_FREE_ALL                             \
  { LG_FREE_WORK }

int LAGraph_MinCut
(
    // outputs
    GrB_Vector* S,
    GrB_Vector* S_bar,
    GrB_Matrix* cut_set,
    // inputs
    LAGraph_Graph G_origin,
    GrB_Matrix R,
    GrB_Index s,
    GrB_Index t,
    char *msg
)
{
  //do a bfs from the source to the sink, stop if the frontier is empty 

  LAGraph_Graph G = NULL;
  GrB_Matrix S_diag=NULL, S_bar_diag=NULL;
  GrB_Index n = 0;
  GrB_Matrix_nrows(&n, R);

  LG_TRY(LAGraph_CheckGraph(G_origin, msg));
  LG_ASSERT (S != NULL, GrB_NULL_POINTER) ;
  LG_ASSERT (S_bar != NULL, GrB_NULL_POINTER) ;
  LG_ASSERT (cut_set != NULL, GrB_NULL_POINTER) ;
  LG_ASSERT (s >= 0 && t >= 0, GrB_INVALID_VALUE) ;

  //  printf("hit\n");
  GrB_Matrix A = G_origin->A;
  
  LG_TRY(GrB_Matrix_new(cut_set, GrB_FP64, n, n));
  LG_TRY(GrB_Vector_new(S_bar, GrB_INT64, n));
  //S is allocated during the bfs
  
  LG_TRY(LAGraph_New(&G, &R, LAGraph_ADJACENCY_DIRECTED, msg));
  LG_TRY(LAGr_BreadthFirstSearch(S, NULL, G, s, msg));


  LG_TRY(GrB_assign(*S_bar, *S, NULL, 1, GrB_ALL, n, GrB_DESC_SC));
  LG_TRY(GrB_assign(*S, *S, NULL, 1, GrB_ALL, n, GrB_DESC_S));

  GRB_TRY(GrB_Matrix_diag(&S_diag, *S, 0));
  GRB_TRY(GrB_Matrix_diag(&S_bar_diag, *S_bar, 0));

  
  GRB_TRY(GrB_mxm(*cut_set, NULL, NULL, GxB_PLUS_TIMES_INT64, A, S_bar_diag, NULL));
  GRB_TRY(GrB_mxm(*cut_set, NULL, NULL, GxB_PLUS_TIMES_INT64, S_diag, *cut_set, NULL));

  
  LG_FREE_ALL;
  return (GrB_SUCCESS);
}
