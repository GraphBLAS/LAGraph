#include <LAGraph.h>
#include "LG_internal.h"
#include <LAGraph.h>

int LAGraph_MinCut
(
    // outputs
    GrB_Vector* S,
    GrB_Vector* S_bar,
    // inputs
    GrB_Matrix R,
    GrB_Index s,
    GrB_Index t,
    char *msg
)
{
  //do a bfs from the source to the sink, stop if the frontier is empty

  LAGraph_Graph G = NULL;
  GrB_Index n = 0;
  GrB_Matrix_nrows(&n, R);
  
  LG_TRY(LAGraph_New(&G, &R, LAGraph_ADJACENCY_DIRECTED, msg));
  LG_TRY(LAGr_BreadthFirstSearch(S, NULL, G, s, msg));

  LG_TRY(GrB_assign(*S_bar, *S, NULL, 1, GrB_ALL, n, GrB_DESC_SC));
  LG_TRY(GrB_assign(*S, *S, NULL, 1, GrB_ALL, n, GrB_DESC_S));

  LAGraph_Delete(&G, msg);
  return (GrB_SUCCESS);
}
