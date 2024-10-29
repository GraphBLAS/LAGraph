#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_Identity
(
    GrB_Matrix I, // //output Identity Matrix
    char* msg
){
    GrB_Index n;
    GRB_TRY(GrB_Matrix_nrows(&n, I));
    for(int i =0;i<n;i++){
        GRB_TRY(GrB_Matrix_setElement_FP64(I,1.0,i,i));
    }
    return 0;
}