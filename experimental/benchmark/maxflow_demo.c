
#include "LAGraphX.h"
#include "../../src/benchmark/LAGraph_demo.h"

#undef LAGRAPH_CATCH

#define LAGRAPH_CATCH(info)                                                    \
  {                                                                            \
    GrB_free(&R);                                                              \
    LAGraph_Delete(&G, msg);                                                   \
    return (info);                                                             \
  }

#define LG_FREE_ALL                                                            \
  { GrB_free(&R); }



int main (int argc, char ** argv){

  char msg[LAGRAPH_MSG_LEN];
  LAGraph_Graph G = NULL;
  GrB_Matrix R = NULL;
  double flow = 0;
  GrB_Index T=0, S=0;

  LAGRAPH_TRY(LAGraph_Init(msg));
  
  //read in graph
  LAGRAPH_TRY(LAGraph_MMRead(&R, stdin, msg));
  GRB_TRY(GrB_Matrix_nrows(&T, R));
  T--;
  LAGRAPH_TRY(readproblem(&G, NULL, false, true, false, NULL, true, argc, argv));
  double time = LAGraph_WallClockTime();
  LAGRAPH_TRY(LAGraph_MaxFlow(G, S, T, &flow, msg));
  time = LAGraph_WallClockTime() - time;

  printf("Time for LAGraph_MaxFlow: %g sec\n", time);
  GrB_free(&R);
  LAGraph_Delete(&G, msg);
  LAGRAPH_TRY(LAGraph_Finalize(msg));

  return GrB_SUCCESS;
  
}
