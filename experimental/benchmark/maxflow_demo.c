
#include "LAGraphX.h"
#include "../../src/benchmark/LAGraph_demo.h"

#undef LAGRAPH_CATCH
#undef GRB_CATCH

#define LAGRAPH_CATCH(info)                                                    \
  {                                                                            \
    LAGraph_Delete(&G, msg);                                                   \
    return (info);                                                             \
  }

#define GRB_CATCH(info) LAGRAPH_CATCH(info)



int main (int argc, char ** argv){

  char msg[LAGRAPH_MSG_LEN];
  LAGraph_Graph G = NULL;

  double flow = 0;
  GrB_Index T=0, S=0;

  LAGRAPH_TRY(LAGraph_Init(msg));

  //read in graph
  LAGRAPH_TRY(readproblem(&G, NULL, false, true, false, NULL, true, argc, argv));
  printf ("read in the problem\n") ;
  double t1 = LAGraph_WallClockTime();
  LAGRAPH_TRY(LAGraph_Cached_AT(G, msg));
  t1 = LAGraph_WallClockTime() - t1 ;
  printf ("cached AT, time %g\n", t1) ;
  t1 = LAGraph_WallClockTime();
  LAGRAPH_TRY(LAGraph_Cached_EMin(G, msg));
  t1 = LAGraph_WallClockTime() - t1 ;
  printf ("cached EMin, time %g\n", t1) ;

  char* end1, *end2;
  S = strtoul(argv[2], &end1, 10);
  T = strtoul(argv[3], &end2, 10);
  if(argc > 4){
    int num_threads = atoi(argv[4]);
    LAGRAPH_TRY(LAGraph_SetNumThreads(1, num_threads, msg));
  }

  int nthreads_outer, nthreads_inner ;
  LAGRAPH_TRY(LAGraph_GetNumThreads(&nthreads_outer, &nthreads_inner, msg));
  printf ("nthreads: %d %d\n", nthreads_outer, nthreads_inner) ;

  if(end1 == 0 || end2 == 0){
    printf("values for source and sink are incorrect.\n");
  }
  printf("Starting max flow from %ld to %ld\n", S, T);

  // LG_SET_BURBLE(1);
  double time = LAGraph_WallClockTime();
  LAGRAPH_TRY(LAGr_MaxFlow(&flow, NULL, G, S, T, msg));
  time = LAGraph_WallClockTime() - time;
  printf("Time for LAGraph_MaxFlow: %g sec\n", time);
  printf("Max Flow is: %lf\n", flow);

  LAGraph_Delete(&G, msg);
  LAGRAPH_TRY(LAGraph_Finalize(msg));

  return GrB_SUCCESS;
}
