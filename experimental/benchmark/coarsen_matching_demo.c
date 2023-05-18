#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

int main(int argc, char **argv)
{
    char msg[1024] ;
    LAGraph_Init (msg) ;
    LAGraph_Random_Init (msg) ;
    GrB_Matrix test = NULL , test2 = NULL ;
    GRB_TRY (LAGraph_Random_Matrix (&test, GrB_BOOL, 3, 5, 0.5, 42, msg)) ;
    GRB_TRY (LAGraph_Random_Matrix (&test2, GrB_BOOL, 5, 3, 0.2, 93, msg)) ;
    GRB_TRY (GrB_transpose (test2, NULL, NULL, test, NULL)) ;
    return (GrB_SUCCESS) ;
}