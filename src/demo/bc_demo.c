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

#include "LAGraph_Test.h"

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 4
// #define THREAD_LIST 40, 20, 10, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&c2) ;                \
    GrB_free (&centrality) ;        \
    GrB_free (&SourceNodes) ;       \
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
    GrB_Vector centrality = NULL, c2 = NULL ;
    GrB_Matrix SourceNodes = NULL ;

    // start GraphBLAS and LAGraph
    LAGraph_TRY (LAGraph_Init (msg)) ;
    bool burble = false ;
    GrB_TRY (GxB_set (GxB_BURBLE, burble)) ;

    int batch_size = 4 ;

    //--------------------------------------------------------------------------
    // determine # of threads to use (TODO: make this a utility
    //--------------------------------------------------------------------------

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
    double tt2 [nthreads_max+1] ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGraph_TRY (LAGraph_Test_ReadProblem (&G, &SourceNodes,
        false, false, true, NULL, false, argc, argv, msg)) ;
    GrB_Index n, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

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
            double tic [2] ;
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

        // TODO: check results

        // GrB_TRY (GxB_print (centrality, 2)) ;
        GrB_free (&centrality) ;

        // if burble is on, just do the first batch
        if (burble) break ;
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
