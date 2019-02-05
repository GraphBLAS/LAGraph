//------------------------------------------------------------------------------
// bfs_test: read in (or create) a matrix and test BFS
//------------------------------------------------------------------------------

// LAGraph, (... list all authors here) (c) 2019, All Rights Reserved.
// http://graphblas.org  See LAGraph/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "LAGraph.h"

#define OK(method)                                                          \
{                                                                           \
    GrB_Info info = method ;                                                \
    if (! (info == GrB_SUCCESS || info == GrB_NO_VALUE))                    \
    {                                                                       \
        fprintf (stderr, "bfs_test failure: [%d] %s\n", info, GrB_error ( )) ;\
        FREE_ALL ;                                                          \
        return (info) ;                                                     \
    }                                                                       \
}

#define FREE_ALL                    \
{                                   \
    GrB_free (&AT) ;                \
    GrB_free (&A) ;                 \
    GrB_free (&v) ;                 \
}

GrB_Matrix AT = NULL ;
GrB_Matrix A = NULL ;
GrB_Matrix Abool = NULL ;
GrB_Vector v = NULL ;

//------------------------------------------------------------------------------
// bfs_get: read in (or create) a matrix
//------------------------------------------------------------------------------

GrB_Info bfs_get
(
    FILE *f,
    GrB_Index n,
    GrB_Index nvals,
    bool make_symmetric,
    bool no_diagonal,
    uint64_t *seed
)
{

    //--------------------------------------------------------------------------
    // read in a matrix from a file, or create a random matrix
    //--------------------------------------------------------------------------

    if (f == NULL)
    {
        // pattern-only
        OK (LAGraph_random (&A, GrB_BOOL, n, n, nvals,
            true, make_symmetric, false, false, no_diagonal, seed)) ;
    }
    else
    {
        // read in the file
        OK (LAGraph_mmread (&A, f)) ;
    }

    // finish any pending computations
    GrB_Matrix_nvals (&nvals, A) ;
    // GxB_fprint (A, GxB_COMPLETE, stdout) ;
    return (GrB_SUCCESS) ;
}



//------------------------------------------------------------------------------
// bfs_test main: test a single matrix using multple methods
//------------------------------------------------------------------------------

// usage:
// bfs_test s < in > out
// s is the staring node, in is the Matrix Market file, out is the level set.

int main (int argc, char **argv)
{

    LAGraph_init ( ) ;

    // read in a matrix from stdin (Matrix Market format)
    OK (LAGraph_mmread (&A, stdin)) ;

    // get the size of the problem.
    GrB_Index nrows, ncols ;
    OK (GrB_Matrix_nrows (&nrows, A)) ;
    OK (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = nrows ;

    // typecast to boolean and change all values to true.
    // TODO: create a Utility function LAGraph_to_pattern that converts
    // a matrix (of any selected type)?
    OK (GrB_Matrix_new (&Abool, GrB_BOOL, nrows, ncols)) ;
    OK (GrB_apply (Abool, NULL, NULL, LAGraph_TRUE_BOOL, A, NULL)) ;
    GrB_free (&A) ;
    A = Abool ;
    Abool = NULL ;

    // GxB_fprint (A, GxB_COMPLETE, stderr) ;

    // get the source node
    GrB_Index s = 0 ;
    if (argc > 1)
    {
        sscanf (argv [1], "%" PRIu64, &s) ; 
    }

    // run the BFS on node s
    GrB_Vector v ;
    OK (LAGraph_bfs_simple (&v, A, s, n)) ;

    // now the BFS on node s using push-pull instead
    GrB_Vector v2 ;
    OK (GrB_Matrix_new (&AT, GrB_BOOL, ncols, nrows)) ;
    OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
    OK (LAGraph_bfs_pushpull (&v2, A, AT, s, n)) ;

    // TODO: you can't typecast from vectors to matrices ...
    // (except in SuiteSparse:GraphBLAS).  Need a function
    // LAGraph_Vector_isequal.
    bool isequal = false ;
    OK (LAGraph_isequal (&isequal, (GrB_Matrix) v, (GrB_Matrix) v2, NULL)) ;

    if (isequal)
    {
        fprintf (stderr, "results of simple and push-pull are OK\n") ;
    }
    else
    {
        fprintf (stderr, "ERROR! simple and push-pull differ\n") ;
        abort ( ) ;
    }

    // write the result to stdout (check them outside of this main program)
    for (int64_t i = 0 ; i < n ; i++)
    {
        int64_t x = 0 ;
        OK (GrB_Vector_extractElement (&x, v, i)) ;
        printf ("%" PRIu64 "\n", x) ;
    }

    FREE_ALL ;
    LAGraph_finalize ( ) ;
    return (GrB_SUCCESS) ;
}

