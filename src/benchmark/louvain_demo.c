#include "LAGraph_demo.h"
#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0
#define LG_FREE_ALL                             \
{                                               \
    GrB_free (&A) ;                             \
    GrB_free (&S) ;                             \
    LAGraph_Delete (&G, msg) ;                  \
}
int main(int argc,char** argv){
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;

    GrB_Matrix A = NULL ;
    GrB_Matrix S = NULL;
    bool burble = false ;
    demo_init (burble) ;
    
    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max, nthreads_outer, nthreads_inner ;
    LAGRAPH_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, msg)) ;
    nthreads_max = nthreads_outer * nthreads_inner ;
    if (Nthreads [1] == 0)
    {
        Nthreads [1] = nthreads_max ;
        for (int t = 2 ; t <= nt ; t++) {
            Nthreads [t] = Nthreads [t-1] / 2 ;
            if (Nthreads [t] == 0) nt = t-1 ;
        }
    }
    printf ("threads to test: ") ;
    for(int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
        printf (" %d", nthreads) ;
    }
    printf ("\n") ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        false, false, true, NULL, false, argc, argv)) ;
    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;
    GrB_Matrix S = NULL; // change to pointer
    double t1 = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads_max, msg)) ;
    LAGRAPH_TRY (LAGraph_Louvain (S,G,msg)) ;
    t1 = LAGraph_WallClockTime ( ) - t1 ;
    printf ("warmup: %10.4f sec\n", t1) ;

    int ntrials = 16 ;
    printf ("# of trials: %d\n", ntrials) ;
    int iters = 0, itermax = 100 ;
    for (int kk = 1 ; kk <= nt ; kk++)
    {
        int nthreads = Nthreads [kk] ;
        if (nthreads > nthreads_max) continue ;
        LAGRAPH_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;
        printf ("\n--------------------------- nthreads: %2d\n", nthreads) ;

        double total_time = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            GrB_free (&S) ;
            double t1 = LAGraph_WallClockTime ( ) ;
            LAGRAPH_TRY (LAGraph_Louvain (S,G,msg)) ;
            t1 = LAGraph_WallClockTime ( ) - t1 ;
            printf ("trial: %2d time: %10.4f sec\n", trial, t1) ;
            total_time += t1 ;
        }

        double t = total_time / ntrials ;
        printf ("GAP: %3d: avg time: %10.3f (sec), "
                "rate: %10.3f iters: %d\n", nthreads,
                t, 1e-6*((double) nvals) * iters / t, iters) ;
        fprintf (stderr, "GAP: Avg: PR %3d: %10.3f sec: %s\n",
             nthreads, t, matrix_name) ;

    }
    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}