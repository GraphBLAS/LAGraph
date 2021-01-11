//------------------------------------------------------------------------------
// test_bc: betweenness centrality for the GAP benchmark
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University

// usage:
// test_bc < matrixfile.mtx
// test_bc matrixfile.mtx sourcenodes.mtx

#include "LAGraph2.h"

#define NTHREAD_LIST 2
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 4
// #define THREAD_LIST 40, 20, 10, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A) ;                 \
    GrB_free (&AT) ;                \
    GrB_free (&Abool) ;             \
    GrB_free (&centrality) ;        \
    GrB_free (&SourceNodes) ;       \
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
    GrB_Matrix AT = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector centrality = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    LAGraph_TRY (LAGraph_Init (msg)) ;
    GrB_TRY (GxB_set (GxB_BURBLE, false)) ;

    uint64_t seed = 1;
    FILE *f ;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max ;
    LAGraph_TRY (LAGraph_GetNumThreads (&nthreads_max, NULL)) ;
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

    double tt [nthreads_max+1] ;

    int batch_size = 4 ;

    //--------------------------------------------------------------------------
    // read in a matrix from a file (TODO: make this a Test2/Utility)
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ; 
    double tic [2] ;
    LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;

    if (argc > 1)
    {
        // Usage:
        //      ./test_bc matrixfile.mtx sources.mtx
        //      ./test_bc matrixfile.grb sources.mtx

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
        }

        // read in source nodes in Matrix Market format from the input file
        if (argc > 2)
        {
            filename = argv [2] ;
            printf ("sources: %s\n", filename) ;
            f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Source node file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGraph_TRY (LAGraph_MMRead (&SourceNodes, f, msg)) ;
            fclose (f) ;
        }
    }
    else
    {

        // Usage:  ./test_bc < matrixfile.mtx
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
        if (G->A_pattern_is_symmetric)
        {
            // if G->A has a symmetric pattern, declare the graph undirected
            // and free G->AT since it isn't needed.  The BFS only looks at
            // the pattern of A anyway.
            G->kind = LAGRAPH_ADJACENCY_UNDIRECTED ;
            GrB_TRY (GrB_Matrix_free (&(G->AT))) ;
        }
    }

    LAGraph_TRY (LAGraph_DisplayGraph (G, 0, msg)) ;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    #define NSOURCES 32

    if (SourceNodes == NULL)
    {
        GrB_TRY (GrB_Matrix_new (&SourceNodes, GrB_INT64, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            int64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
            // SourceNodes [k] = i 
            LAGraph_TRY (GrB_Matrix_setElement (SourceNodes, i, k, 0)) ;
        }
    }

    GrB_Index nsource ;
    GrB_TRY (GrB_Matrix_nrows (&nsource, SourceNodes)) ;
    if (nsource % batch_size != 0)
    {
        printf ("SourceNode size must be multiple of batch_size (%d)\n",
            batch_size) ;
        exit (1) ;
    }

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    int ntrials = 0 ;

    for (int64_t kstart = 0 ; kstart < nsource ; kstart += batch_size)
    {

        //----------------------------------------------------------------------
        // Create batch of vertices to use in traversal
        //----------------------------------------------------------------------

        ntrials++ ;
        printf ("\nTrial %d : sources: [", ntrials) ;
        GrB_Index vertex_list [batch_size] ;
        for (int64_t k = 0 ; k < batch_size ; k++)
        {
            // get the kth source node
            int64_t source = -1 ;
            GrB_TRY (GrB_Matrix_extractElement (&source, SourceNodes,
                k + kstart, 0)) ;
            // subtract one to convert from 1-based to 0-based
            source-- ;
            vertex_list [k] = source  ;
            printf (" %"PRIu64, source) ;
        }
        printf (" ]\n") ;

        //----------------------------------------------------------------------
        // Compute betweenness centrality using batch algorithm
        //----------------------------------------------------------------------

        // back to default
        GxB_set (GxB_NTHREADS, nthreads_max) ;

        for (int t = 1 ; t <= nt ; t++)
        {
            if (Nthreads [t] > nthreads_max) continue ;
            GrB_TRY (GxB_set (GxB_NTHREADS, Nthreads [t])) ;
            GrB_free (&centrality) ;
            LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
            LAGraph_TRY (LAGraph_VertexCentrality_Betweenness
                (&centrality, G, vertex_list, batch_size, msg)) ;

            double t2 ;
            LAGraph_TRY (LAGraph_Toc (&t2, tic, msg)) ;
            printf ("BC time %2d: %12.4f (sec)\n", Nthreads [t], t2) ;
            fflush (stdout) ;
            tt [t] += t2 ;
        }

        //----------------------------------------------------------------------
        // check result
        //----------------------------------------------------------------------

        // TODO

        // GrB_TRY (GxB_print (centrality, 2)) ;
        GrB_free (&centrality) ;

        // HACK: uncomment this to just do the first batch
        // break ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    printf ("\nntrials: %d\n", ntrials) ;

    printf ("\n") ;
    for (int t = 1 ; t <= nt ; t++)
    {
        if (Nthreads [t] > nthreads_max) continue ;
        double t2 = tt [t] / ntrials ;
        printf ("Ave BC %2d: %10.3f sec, rate %10.3f\n",
            Nthreads [t], t2, 1e-6*((double) nvals) / t2) ;
        fprintf (stderr, "Avg: BC %3d: %10.3f sec: %s\n",
            Nthreads [t], t2, matrix_name) ;
    }

    LAGRAPH_FREE_ALL;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}

