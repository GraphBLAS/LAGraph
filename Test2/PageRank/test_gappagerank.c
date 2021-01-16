//------------------------------------------------------------------------------
// test_gappagerank: read in (or create) a matrix and test the GAP PageRank
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M and Gabor Szarnyas, BME

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#include "LAGraph2.h"

#define LAGRAPH_FREE_ALL                        \
{                                               \
    GrB_free (&A) ;                             \
    GrB_free (&Abool) ;                         \
    GrB_free (&PR) ;                            \
    LAGraph_Delete (&G, msg) ;                  \
    if (f != NULL) fclose (f) ;                 \
}

#define LAGraph_CATCH(status)                                               \
{                                                                           \
    printf ("LAGraph error: %s line: %d, status: %d: %s\n", __FILE__,       \
        __LINE__, status, msg) ;                                            \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}

#define GrB_CATCH(info)                                                     \
{                                                                           \
    printf ("GraphBLAS error: %s line: %d, info: %d: %s\n", __FILE__,       \
        __LINE__, info, msg) ;                                              \
    LAGRAPH_FREE_ALL ;                                                      \
    return (-1) ;                                                           \
}

int main (int argc, char **argv)
{

    printf ("%s v%d.%d.%d [%s]\n",
        GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR,
        GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB,
        GxB_IMPLEMENTATION_DATE) ;

    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;

    GrB_Matrix A = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector PR = NULL ;
    FILE *f = NULL ;

    // start GraphBLAS and LAGraph
    LAGraph_TRY (LAGraph_Init (msg)) ;
    GrB_TRY (GxB_set (GxB_BURBLE, false)) ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ;
    GrB_TRY (GxB_get (GxB_NTHREADS, &nthreads_max)) ;
    if (Nthreads [1] == 0)
    {
        // create thread list automatically
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++)
        {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test: ") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file (TODO: make this a Test2/Utility)
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ; 
    double tic [2] ;
    LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;

    if (argc > 1)
    {
        // Usage:
        //      ./test_gappagerank matrixfile.mtx sources.mtx
        //      ./test_gappagerank matrixfile.grb sources.mtx

        // read in the file in Matrix Market format from the input file
        char *filename = argv [1] ;
        printf ("matrix: %s\n", filename) ;

        // find the filename extension
        size_t len = strlen (filename) ;
        char *ext = NULL ;
        for (int k = len-1 ; k >= 0 ; k--)
        {
            if (filename [k] == '.')
            {
                ext = filename + k ;
                printf ("[%s]\n", ext) ;
                break ;
            }
        }
        bool is_binary = (ext != NULL && strncmp (ext, ".grb", 4) == 0) ;

        if (is_binary)
        {
            printf ("Reading binary file: %s\n", filename) ;
            LAGraph_TRY (LAGraph_BinRead (&A, filename, msg)) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGraph_TRY (LAGraph_MMRead (&A, f, msg)) ;
            fclose (f) ;
            f = NULL ;
        }

    }
    else
    {

        // Usage:  ./test_gappagerank < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGraph_TRY (LAGraph_MMRead (&A, stdin, msg)) ;
    }

    //--------------------------------------------------------------------------
    // convert to boolean, pattern-only
    //--------------------------------------------------------------------------

    LAGraph_TRY (LAGraph_Pattern (&Abool, A, msg)) ;
    GrB_free (&A) ;
    A = Abool ;
    Abool = NULL ;
    GrB_TRY (GrB_wait (&A)) ;

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    GrB_Index n = nrows ;
    if (nrows != ncols) { printf ("A must be square\n") ; abort ( ) ; }
    double t_read ;
    LAGraph_TRY (LAGraph_Toc (&t_read, tic, msg)) ;
    printf ("read time: %g\n", t_read) ;

    //--------------------------------------------------------------------------
    // construct the graph
    //--------------------------------------------------------------------------

    bool A_is_symmetric =
        (nrows == 134217726 ||  // HACK for kron
         nrows == 134217728) ;  // HACK for urand

    if (A_is_symmetric)
    {
        // A is known to be symmetric
        // TODO: LAGraph_New should set G->A_pattern_is_symmetric if
        // the G->kind is LAGRAPH_ADJACENCY_UNDIRECTED
        LAGraph_TRY (LAGraph_New (&G, &A, LAGRAPH_ADJACENCY_UNDIRECTED, false,
            msg)) ;
        G->A_pattern_is_symmetric = true ;
    }
    else
    {
        // compute G->AT and determine if A has a symmetric pattern
        LAGraph_TRY (LAGraph_New (&G, &A, LAGRAPH_ADJACENCY_DIRECTED, false,
            msg)) ;
        LAGraph_TRY (LAGraph_Property_ASymmetricPattern (G, msg)) ;
        A_is_symmetric = G->A_pattern_is_symmetric ;
        if (A_is_symmetric)
        {
            // if G->A has a symmetric pattern, declare the graph undirected
            // and free G->AT since it isn't needed.  The BFS only looks at
            // the pattern of A anyway.
            G->kind = LAGRAPH_ADJACENCY_UNDIRECTED ;
            GrB_TRY (GrB_Matrix_free (&(G->AT))) ;
        }
    }

    LAGraph_TRY (LAGraph_Property_RowDegree (G, msg)) ;

    //--------------------------------------------------------------------------
    // compute the pagerank
    //--------------------------------------------------------------------------

    // the GAP benchmark requires 16 trials
    int ntrials = 16 ;
    // ntrials = 1 ;    // HACK to run just one trial
    printf ("# of trials: %d\n", ntrials) ;

    float damping = 0.85 ;
    float tol = 1e-4 ;
    int iters = 0, itermax = 100 ;

    for (int kk = 1 ; kk <= nt ; kk++)
    {
        int nthreads = Nthreads [kk] ;
        if (nthreads > nthreads_max) continue ;
        GxB_set (GxB_NTHREADS, nthreads) ;
        printf ("\n--------------------------- nthreads: %2d\n", nthreads) ;

        double total_time = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            GrB_free (&PR) ;
            LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
            LAGraph_TRY (LAGraph_VertexCentrality_PageRankGAP (&PR, G,
                damping, tol, itermax, &iters, msg)) ;
            double t1 ;
            LAGraph_TRY (LAGraph_Toc (&t1, tic, NULL)) ;
            printf ("trial: %2d time: %10.4f sec\n", trial, t1) ;
            total_time += t1 ;
        }

        double t = total_time / ntrials ;
        printf ("3f:%3d: avg time: %10.3f (sec), "
                "rate: %10.3f iters: %d\n", nthreads,
                t, 1e-6*((double) nvals) * iters / t, iters) ;
        fprintf (stderr, "Avg: PR (3f)      %3d: %10.3f sec: %s\n",
             nthreads, t, matrix_name) ;

    }
    
    //--------------------------------------------------------------------------
    // check result
    //--------------------------------------------------------------------------

    // GxB_print (PR, 2) ;
    // TODO

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_ALL ;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}

