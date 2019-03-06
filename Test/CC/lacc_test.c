#include "LAGraph.h"

#define LAGRAPH_FREE_ALL    \
{                           \
    GrB_free (&A) ;         \
}


int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL ;
    GrB_init (GrB_NONBLOCKING) ;
    // self edges are OK

    FILE *f ;
    if (argc == 1)
    {
        f = stdin ; 
    }
    else
    {
        f = fopen (argv[0], "r") ;
        if (f == NULL)
        {
            printf ("unable to open file [%s]\n", argv[0]) ;
            return (GrB_INVALID_VALUE) ;
        }
    }

    LAGRAPH_OK (LAGraph_mmread (&A, f)) ;

    LAGRAPH_OK (LAGraph_lacc (A)) ;

    LAGRAPH_FREE_ALL ;
    LAGRAPH_OK (GrB_finalize ( )) ;
}

