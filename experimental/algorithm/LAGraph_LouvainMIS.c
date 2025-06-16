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
        GrB_set (GrB_GLOBAL,true, GxB_BURBLE);
    // Shortened monoids, Binary ops, and Semirings
    GrB_Monoid plusmon = GrB_PLUS_MONOID_FP64;
    GrB_Monoid timesmon = GrB_TIMES_FP64;
    GrB_BinaryOp plusf64 = GrB_PLUS_FP64;
    GrB_BinaryOp timesf64 = GrB_TIMES_FP64;
    GrB_BinaryOp divf64 = GrB_DIV_FP64;
    GrB_Semiring stdmxm = GrB_PLUS_TIMES_SEMIRING_FP64;

    // GrB_BinaryOp accum_S_vector = NULL;
    // Declarations
    GrB_Vector iset = NULL;
    GrB_Vector k = NULL;
    GrB_Vector x = NULL;
    GrB_Vector y = NULL;
    GrB_Vector neighbours = NULL;
    GrB_Vector S_vector = NULL;
    GrB_Matrix S = NULL;
    GxB_Container S_container = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix A_iset = NULL;
    GrB_Index n ; 
    GrB_Matrix W;
    // -------------------------------------------------------//

    // -------------------------------------------------------//
    
    // Initializing 
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT(S_result != NULL, GrB_NULL_POINTER);
    
    A = G->A;
    dbg(A);
    
    //k = [+_j A(:,j)]
    GRB_TRY(GrB_Matrix_nrows(&n,A));
    GRB_TRY(GrB_Vector_new(&k, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&y, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&x,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&neighbours,GrB_FP64,n));
    GRB_TRY(GrB_Vector_new(&S_vector,GrB_FP64,n));
    GRB_TRY(GrB_Matrix_new(&W,GrB_FP64,n,n));
    GRB_TRY(GrB_Matrix_new(&A_iset,GrB_FP64,n,n));
    GRB_TRY(GxB_Container_new(&S_container));


    GRB_TRY (GrB_Matrix_reduce_Monoid(k, NULL, NULL,plusmon, A, NULL));
    // dbg(k);
    double m;
    GRB_TRY(GrB_Vector_reduce_FP64(&m,NULL,plusmon,k,NULL));
    m*=.5;
    printf("%ld",m);
    
    // S <- I
    GRB_TRY(GrB_assign (x, NULL, NULL, 1, GrB_ALL, n, NULL)) ;
    // GxB_print(i,5);

    GRB_TRY(GrB_Matrix_diag(&S,x,0));
    GRB_TRY(GrB_set(S,false,GxB_ISO));
    GrB_set (S, GxB_SPARSE, GxB_SPARSITY_CONTROL) ;
    // dbg(S);
    GRB_TRY(GxB_unload_Matrix_into_Container(S,S_container,NULL));
    // dbg(S_container->x);
    GRB_TRY(GrB_Vector_setElement(S_container->i,33,31));
    GRB_TRY(GrB_Vector_setElement(S_container->i,32,30));
    // GRB_TRY(GrB_Vector_setElement_BOOL(S_container->x,false,1));
    GRB_TRY(GxB_load_Matrix_from_Container(S,S_container,NULL));
    dbg(S);
    // Compute Isolate Set
    GRB_TRY(LAGraph_IsolateSets(&iset,G,32121,msg));
    dbg(iset);

    // Compute max change in modularity for each node in the isolate set


    GRB_TRY(GrB_mxm(W,NULL,NULL,GrB_PLUS_TIMES_SEMIRING_FP64,A_iset,S,NULL));
    GRB_TRY(GrB_vxm(y,NULL,NULL,GrB_PLUS_TIMES_SEMIRING_FP64,y,S,GrB_DESC_T0));
    // GRB_TRY(GrB_mxv(,NULL,NULL,GrB_PLUS_TIMES_SEMIRING_FP64,W,y,NULL)); 
    // GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(coef,NULL,NULL,timesf64,k,iset,NULL));
    // GRB_TRY(GrB_Vector_apply_BinaryOp2nd_FP64(coef,NULL,NULL,divf64,coef,m,GrB_DESC_T0));
    // dbg(coef);

    // GRB_TRY(GrB_vxm(neighbours,NULL,NULL,GrB_PLUS_TIMES_SEMIRING_FP64,coef,A,NULL));    //This gives us the neighbours of each node in the Isolate set
    // dbg(neighbours);

    // GRB_TRY(GrB_mxv(c,NULL,NULL,GrB_PLUS_TIMES_SEMIRING_FP64,S,neighbours,NULL)); //This gives all nodes in the S_vectorunity in each of the neighbours S_vectorunity
    // dbg(c);
    // GRB_TRY(GrB_Vector_eWiseMult_BinaryOp(ck,NULL,NULL,timesf64,k,c,NULL));    
    // dbg(ck); 
    


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