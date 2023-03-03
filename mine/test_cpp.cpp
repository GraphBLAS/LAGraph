extern "C" {
   #include "LAGraph.h"
   #include "LAGraphX.h"
   #include "GraphBLAS.h"
}

int main(){
   char msg [1024] ;
   LAGraph_Init (msg) ;
   LAGraph_Random_Init (msg) ;
   GrB_Matrix A ;
   // GrB_Matrix_new (&A, GrB_BOOL, 10, 10) ;
   LAGraph_Random_Matrix (&A, GrB_BOOL, 10, 10, 0.5, 42, msg) ;
   GxB_Matrix_fprint (A, "A", GxB_COMPLETE, stdout) ;
   printf("hi\n");
   GrB_Matrix_free (&A) ;
   LAGraph_Random_Finalize (msg) ;
   LAGraph_Finalize (msg) ;
}
