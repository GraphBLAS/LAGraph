//------------------------------------------------------------------------------
// LAGraph/Test2/ConnectedComponents/test_cc: test LAGraph_ConnectedComponents
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

// Contributed by Tim Davis, Texas A&M

// Usage: test_cc can be used with both stdin or a file as its input,
// in either grb or mtx format.

//------------------------------------------------------------------------------

#include "LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_alg_internal.h"

#undef  LAGraph_FREE_ALL
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
    GrB_Vector components = NULL, components2 = NULL, components3 = NULL,
        components4 = NULL ;
    GrB_Vector components5 = NULL ;

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
    if (readproblem (&G,
        NULL,   // no source nodes
        true,   // make the graph undirected, and symmetrize the matrix
        false,  // do not remove self-edges
        true,   // structural only, no values needed
        NULL,   // no type preference
        false,  // do not ensure all entries positive
        argc, argv) != 0) ERROR ;
    GrB_Index n, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // begin tests
    //--------------------------------------------------------------------------

    double tic [2], t1, t2, t3, t4, t5 ;

    #define NTRIALS 16
    printf ("# of trials: %d\n\n", NTRIALS) ;

    GrB_Index nCC, nCC_first = -1, nCC2 = 0 , nCC3, nCC4 = 0, nCC5 = 0 ;
    for (int trial = 1 ; trial <= nt ; trial++)
    {
        int nthreads = Nthreads [trial] ;
        if (nthreads > nthreads_max) continue ;
        LAGraph_TRY (LAGraph_SetNumThreads (nthreads, NULL)) ;

        t1 = 0 ;
        t2 = 0 ;
        t3 = 0 ;
        t4 = 0 ;
        t5 = 0 ;

        for (int k = 0 ; k < NTRIALS ; k++)
        {

            double ttrial, tcheck ;
            int status ;

            //------------------------------------------------------------------
            // LAGraph_ConnectedComponents
            //------------------------------------------------------------------

            LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
            status = LAGraph_ConnectedComponents (&components, G, msg) ;
            LAGraph_TRY (status) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, NULL)) ;
            t1 += ttrial ;
            printf ("SV6:  trial: %2d time: %10.4f sec\n", k, ttrial) ;
            nCC = countCC (components, n) ;

            //------------------------------------------------------------------
            // LG_CC_FastSV5: slightly faster (using 32-bit integers)
            //------------------------------------------------------------------

            #if LG_SUITESPARSE
            LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
            status = LG_CC_FastSV5 (&components4, G, msg) ;
            LAGraph_TRY (status) ;
            LAGraph_TRY (LAGraph_Toc (&ttrial, tic, NULL)) ;
            t4 += ttrial ;
            printf ("SV5b: trial: %2d time: %10.4f sec\n", k, ttrial) ;
            nCC4 = countCC (components4, n) ;
            #endif

            //------------------------------------------------------------------
            // LG_CC_Boruvka
            //------------------------------------------------------------------

            if (k == 0)
            {
                LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
                status = LG_CC_Boruvka (&components2, G, msg) ;
                LAGraph_TRY (status) ;
                LAGraph_TRY (LAGraph_Toc (&ttrial, tic, NULL)) ;
                t2 += ttrial ;
                printf ("Boru: trial: %2d time: %10.4f sec\n", k, ttrial) ;
                nCC2 = countCC (components2, n) ;
                printf ("nCC2 %ld\n", nCC2) ;
            }

            //------------------------------------------------------------------
            // LAGraph_cc_boruvka
            //------------------------------------------------------------------

            #if 0
            if (k == 0)
            {
                LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
                status = LAGraph_cc_boruvka (&components5, G->A, false) ;
                LAGraph_TRY (status) ;
                LAGraph_TRY (LAGraph_Toc (&ttrial, tic, NULL)) ;
                t5 += ttrial ;
                printf ("oldB: trial: %2d time: %10.4f sec\n", k, ttrial) ;
                nCC5 = countCC (components5, n) ;
                printf ("nCC5 %ld\n", nCC5) ;
            }
            #endif

            //------------------------------------------------------------------
            // LAGraph_cc_lacc
            //------------------------------------------------------------------

            if (k == 0)
            {
                LAGraph_TRY (LAGraph_Tic (tic, NULL)) ;
                status = (LAGraph_cc_lacc (&components3, G->A, false)) ;
                LAGraph_TRY (status) ;
                LAGraph_TRY (LAGraph_Toc (&ttrial, tic, NULL)) ;
                t3 += ttrial ;
                printf ("LACC: trial: %2d time: %10.4f sec\n", k, ttrial) ;
                nCC3 = countCC (components3, n) ;
            }

            //------------------------------------------------------------------
            // check results
            //------------------------------------------------------------------

            if (k == 0)
            {
                // check the result
                printf ("nCC: %g\n", (double) nCC) ;
                nCC_first = nCC ;
                if (nCC2 != nCC) printf ("failure: nCC2 %g err %g\n",
                    (double) nCC2, (double )(nCC - nCC2)) ;
//              if (nCC5 != nCC) printf ("failure: nCC5 %g err %g\n",
//                  (double) nCC5, (double )(nCC - nCC5)) ;
                if (nCC3 != nCC) printf ("failure: nCC3 %g\n", (double) nCC3) ;
                if (nCC4 != nCC) printf ("failure: nCC4 %g\n", (double) nCC4) ;
#if 0 & LG_CHECK_RESULT
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
                    printf ("test failure: # components differs %g %g\n",
                        (double) nCC, (double) nCC_first) ;
                }
            }
            GrB_free (&components) ;
            GrB_free (&components2) ;
            GrB_free (&components3) ;
            GrB_free (&components4) ;
            GrB_free (&components5) ;
        }

        double ttt = t1 / NTRIALS ;
        printf("LG2:SV6:  threads: %2d time: %10.4f  # of CC: %g\n\n",
            nthreads, ttt, (double) nCC) ;
        fprintf (stderr,
            "Avg: CC (sv6)     %3d: %10.3f sec: %s\n",
            nthreads, ttt, matrix_name) ;

        ttt = t4 / NTRIALS ;
        printf("LG2:SV5b: threads: %2d time: %10.4f  # of CC: %g\n\n",
            nthreads, ttt, (double) nCC4) ;
        fprintf (stderr,
            "Avg: CC (sv5b)    %3d: %10.3f sec: %s\n",
            nthreads, ttt, matrix_name) ;

        ttt = t2 ;
        if (ttt > 0)
        {
            printf("LG2:Boru: threads: %2d time: %10.4f  # of CC: %g\n\n",
                nthreads, ttt, (double) nCC2) ;
            fprintf (stderr,
                "Avg: CC (Boruvka) %3d: %10.3f sec: %s\n",
                nthreads, ttt, matrix_name) ;
        }

        ttt = t5 ;
        if (ttt > 0)
        {
            printf("old:boru: threads: %2d time: %10.4f  # of CC: %g\n\n",
                nthreads, ttt, (double) nCC5) ;
            fprintf (stderr,
                "Avg: CC old_boru  %3d: %10.3f sec: %s\n",
                nthreads, ttt, matrix_name) ;
        }

        ttt = t3 ;
        if (ttt > 0)
        {
            printf("LG2:LACC: threads: %2d time: %10.4f  # of CC: %g\n\n",
                nthreads, ttt, (double) nCC3) ;
            fprintf (stderr,
                "Avg: CC (LACC)    %3d: %10.3f sec: %s\n",
                nthreads, ttt, matrix_name) ;
        }
        printf("\n");
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    LAGraph_FREE_ALL;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}
