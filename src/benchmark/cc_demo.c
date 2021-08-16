//------------------------------------------------------------------------------
// LAGraph/Test2/ConnectedComponents/test_cc: test LAGraph_ConnectedComponents
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage: test_cc can be used with both stdin or a file as its input,
// in either grb or mtx format.

#include "LAGraph_demo.h"

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&components) ;        \
}

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

GrB_Index countCC (GrB_Vector f, GrB_Index n)
{
    GrB_Index nCC = 0;
    GrB_Index *w_val = (GrB_Index *) LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    if (w_val == NULL) { printf ("out of memory\n") ; abort ( ) ; }
    // FIXME: NULL parameter to GrB_Vector_extractTuples is SS extension.
    GrB_Vector_extractTuples (NULL, w_val, &n, f) ; // FIXME
    for (GrB_Index i = 0; i < n; i++)
    {
        if (w_val[i] == i)
        {
            nCC++ ;
        }
    }
    LAGraph_Free ((void **) &w_val) ;
    return nCC;
}

int main (int argc, char **argv)
{

    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Vector components = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

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

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    if (readproblem (&G, NULL,
        true, false, true, NULL, false, argc, argv) != 0) ERROR ;
    GrB_Index n, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // begin tests
    //--------------------------------------------------------------------------

    double tic [2], t1, t2 ;

    #define NTRIALS 16
    printf ("# of trials: %d\n\n", NTRIALS) ;

    GrB_Index nCC, nCC_first = -1 ;
    for (int trial = 1 ; trial <= nt ; trial++)
    {
        int nthreads = Nthreads [trial] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_TRY (LAGraph_SetNumThreads (nthreads, NULL)) ;

        t1 = 0 ;
        for (int k = 0 ; k < NTRIALS ; k++)
        {
            LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
            int status = (LAGraph_ConnectedComponents (&components, G, msg)) ;
            // printf ("status : %d msg: %s\n", status, msg) ;
            LAGraph_TRY (status) ;
            double ttrial, tcheck ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, NULL)) ;
            t1 += ttrial ;
            printf ("SV5b: trial: %2d time: %10.4f sec\n", k, ttrial) ;
            nCC = countCC (components, n) ;
            if (k == 0)
            {
                // check the result
                printf ("nCC: %ld\n", nCC) ;
                nCC_first = nCC ;
#if LG_CHECK_RESULT
                LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
                int result = LG_check_cc (components, G, msg) ;
                if (result != 0)
                {
                    printf ("test failure: (%d) %s\n", result, msg) ;
                }
                LAGraph_TRY (LAGraph_Toc (&tcheck, tic, NULL)) ;
                LAGraph_TRY (result) ;
                printf ("LG_check_cc passed, time: %g\n", tcheck) ;
#endif
            }
            else
            {
                if (nCC != nCC_first)
                {
                    printf ("test failure: # components differs %ld %ld\n",
                        nCC, nCC_first) ;
                }
            }
            GrB_free (&components) ;
        }
        double ttt = t1 / NTRIALS ;
        printf("LG2:SV5b: threads: %2d time: %10.4f  # of CC: %lu\n\n",
            nthreads, ttt, nCC) ;
        fprintf (stderr,
            "Avg: CC (sv5b.2)  %3d: %10.3f sec: %s\n",
            nthreads, ttt, matrix_name) ;

        printf("\n");
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGraph_FREE_ALL;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}
