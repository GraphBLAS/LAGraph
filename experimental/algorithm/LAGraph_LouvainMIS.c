#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int LAGraph_LouvainMIS(
    //output
    GrB_Matrix *S_result,
    //input 
    LAGraph_Graph G,
    char* msg
){
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;

    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;

    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;
    
}