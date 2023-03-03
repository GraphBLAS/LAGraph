// extern "C" {
   #include "GraphBLAS.h"
   #include "LAGraph.h"
   #include "LAGraphX.h"
// }

int main(){

   GrB_init (GrB_NONBLOCKING) ;
   GxB_Global_Option_set (GxB_BURBLE, true) ;
   char msg [1024] ;
   GrB_Matrix A = NULL ;
// GrB_Matrix_new (&A, GrB_BOOL, 10, 10) ;
   GrB_Info info = LAGraph_Random_Matrix (&A, GrB_FP32, 10, 10, 0.5, 42, msg) ;
   printf ("result %d\n", info) ;
   GxB_Matrix_fprint (A, "A", GxB_COMPLETE, stdout) ;
   GrB_Matrix_free (&A) ;
   GrB_finalize () ;
}

