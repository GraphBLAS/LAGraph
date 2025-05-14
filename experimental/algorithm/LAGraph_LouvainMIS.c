#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define DEBUG 1
#define dbg(x) if (DEBUG) GxB_print(x,5)
#undef LG_FREE_ALL
#define LG_FREE_ALL\
{   \
}
int LAGraph_LouvainMIS(
    //output
    GrB_Matrix *S_result,
    //input 
    LAGraph_Graph G,
    char* msg
){
    LG_CLEAR_MSG ;

    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    if (DEBUG)
        GrB_set (GrB_GLOBAL, true, GxB_BURBLE);
    // Shortened monoids, Binary ops, and Semirings
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;
    GrB_Monoid timesmon = GrB_TIMES_FP64;
    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;
    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;


    // Declarations
    GrB_Vector iset = NULL;
    GrB_Vector v = NULL;
    GrB_Vector k = NULL;
    GrB_Vector x = NULL;
    GrB_Vector coef = NULL;
    GrB_Vector z = NULL;
    GrB_Matrix S = NULL;
    GrB_Matrix A = NULL;
    GrB_Index n ; 

    
    // Initializing 
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT(S_result != NULL, GrB_NULL_POINTER);
    
    A = G->A;
    dbg(A);
    
    //k = [+_j A(:,j)]
    GRB_TRY(GrB_Matrix_nrows(&n,A));
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&coef,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&x,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&z,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&v,GrB_FP64,n));


    double m;
    GRB_TRY(GrB_Vector_reduce_FP64(&m,NULL,plusmon,k,NULL));
    m*=.5;
    printf("%ld",m);
    
    // S <- I
    GRB_TRY(GrB_assign (x, NULL, NULL, 1, GrB_ALL, n, NULL)) ;
    // GxB_print(i,5);
    GRB_TRY(GrB_Matrix_diag(&S,x,0));
    dbg(S);
    // Compute Isolate Set
    GRB_TRY(LAGraph_IsolateSets(&iset,G,123,msg));
    dbg(iset);

    // Compute max change in modularity for each node in the isolate set
    GRB_TRY(GrB_mxv(z,NULL,NULL,stdmxm,S,k,NULL));
    dbg(z);
    GRB_TRY (GrB_Matrix_reduce_Monoid(k, NULL, NULL,plusmon, A, NULL));
    dbg(k);
    GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(coef,NULL,NULL,timesf64,k,iset,NULL));
    
    dbg(coef);

    // Aggregate Graph
    // Iterate till no change is detected
    
    // double Q;
    // double gamma = 1;
    // GRB_TRY(LAGr_Modularity2(&Q,gamma,A,S,msg));
    // // printf("Iterations: %d\n", iter);
    // printf("Q:%.15g\n",Q);
    (*S_result) = S ;
    S = NULL;
    LG_FREE_ALL;
    return 0;

}