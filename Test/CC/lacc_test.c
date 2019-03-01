#include "LAGraph.h"


int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL ;
    GrB_init (GrB_NONBLOCKING) ;
    // self edges are OK

    f = fopen (argv[0], "r") ;
    if (f == NULL)
    {
        printf ("unable to open file [%s]\n", argv[0]) ;
        FREE_ALL ;
        return (GrB_INVALID_VALUE) ;
    }
    OK (LAGraph_mmread (&A, f)) ;
 
    LACC(A);
    
    OK(GrB_finalize ( ) );
}

