
#include "LG_internal.h"
#include <LG_Xtest.h>

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK				\
{						\
 GrB_free(&flow_r);				\
 GrB_free(&flow_c);				\
 GrB_free(&result_vec);				\
}


#define LG_FREE_ALL				\
{						\
  LG_FREE_WORK;					\
}


int LG_check_flow(const GrB_Matrix flow_mtx, char* msg)
{
  GrB_Vector flow_r=NULL, flow_c=NULL, result_vec=NULL ;
  GrB_Index n ;
  double net_flow = -1;
  LG_TRY(GrB_Matrix_nrows(&n, flow_mtx));
  LG_TRY(GrB_Vector_new(&flow_r, GrB_FP64, n));
  LG_TRY(GrB_Vector_new(&flow_c, GrB_FP64, n));
  LG_TRY(GrB_Vector_new(&result_vec, GrB_FP64, n));
  LG_TRY(GrB_reduce(flow_c, NULL, NULL, GrB_PLUS_MONOID_FP64, flow_mtx, NULL));
  LG_TRY(GrB_reduce(flow_r, NULL, NULL, GrB_PLUS_MONOID_FP64, flow_mtx, GrB_DESC_T1));
  LG_TRY(GrB_eWiseAdd(result_vec, NULL, NULL, GrB_MINUS_FP64, flow_r, flow_c, NULL));
  LG_TRY(GrB_reduce(&net_flow, NULL, GrB_PLUS_MONOID_FP64, result_vec, NULL));
  LG_ASSERT_MSG(net_flow == 0, GrB_INVALID_VALUE, "Flow conservation is not followed");
  LG_FREE_WORK ;
  return GrB_SUCCESS;
}
