// extern "C" {
   #include "GraphBLAS.h"
   #include "LAGraph.h"
   #include "LAGraphX.h"
// }

#define OK(method)                                                  \
{                                                                   \
    GrB_Info info = method ;                                        \
    if (!(info == GrB_SUCCESS || info != GrB_NO_VALUE))             \
    {                                                               \
        printf ("error! line %d info %d\n", __LINE__, info) ;       \
        abort ( ) ;                                                 \
    }                                                               \
}

int main(){

   char msg [1024] ;
   GrB_init (GrB_NONBLOCKING) ;
   OK (LAGraph_Random_Init (msg)) ;
   GxB_Global_Option_set (GxB_BURBLE, true) ;
   GrB_Matrix A = NULL ;
// GrB_Matrix_new (&A, GrB_BOOL, 10, 10) ;
   GrB_Info info = LAGraph_Random_Matrix (&A, GrB_FP32, 10, 10, 0.5, 42, msg) ;
   printf ("result %d\n", info) ;
   printf ("message: %s\n", msg) ;
   GxB_Matrix_fprint (A, "A", GxB_COMPLETE, stdout) ;
   GrB_Matrix_free (&A) ;
   OK (LAGraph_Random_Finalize (msg)) ;
   GrB_finalize () ;
}

