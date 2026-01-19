#include <LAGraph.h>
#include "LG_internal.h"
#include <LAGraph.h>

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK				\
{						\
 GrB_free(&frontier);				\
}

#define LG_FREE_ALL				\
  LG_FREE_WORK

int LAGraph_MinCut
(
    // outputs
    GrB_Vector S,
    GrB_Vector S_bar,
    // inputs
    GrB_Matrix R,
    GrB_Index s,
    GrB_Index t,
    char *msg
)
{
  //do a bfs from the source to the sink, stop if the frontier is empty
  //assign the frontier to S
  GrB_Vector frontier = NULL;
  GrB_Index n = 0, n_frontier = 1;

  GRB_TRY(GrB_Matrix_nrows(&n, R));
  GRB_TRY(GrB_Vector_new(&frontier, GrB_INT64, n));
  GRB_TRY(GrB_Vector_setElement(frontier, 1, s));

  //initial assign to S
  GRB_TRY(GrB_assign(S, NULL, NULL, frontier, GrB_ALL, NULL));
  while(n_frontier > 0){
    GRB_TRY(GrB_mxv(frontier, NULL, NULL, GxB_ANY_PAIR_INT64, R, frontier, GrB_DESC_R));
    GRB_TRY(GrB_assign(S, S, NULL, frontier, GrB_ALL, GrB_DESC_SC));
    GRB_TRY(GrB_Vector_nvals(&n_frontier, frontier));
  }

  GRB_TRY(GrB_assign(S_bar, S, NULL, 1, GrB_ALL, n, GrB_DESC_SC));

  LG_FREE_ALL;
  return (GrB_SUCCESS);
}
