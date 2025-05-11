#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

#define VERBOSE

#define DEFAULT_SIZE 100
#define DEFAULT_DENSITY 0.5
#define DEFAULT_SEED 42

#define NTHREAD_LIST 0
#define THREAD_LIST 8

// #define SHOW_RESULTS

#define LG_FREE_ALL                     \
{                                       \
    LAGraph_Delete (&G, NULL) ;         \
    GrB_free (&A) ;                     \
    GrB_free (&coarsened) ;             \
    GrB_free (&parent_result) ;         \
    GrB_free (&newlabel_result) ;       \
    GrB_free (&inv_newlabel_result) ;   \
}                                       \

int main(int argc, char **argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Matrix coarsened = NULL ;
    GrB_Vector parent_result = NULL, newlabel_result = NULL, inv_newlabel_result = NULL ;

    bool burble = false ; 
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in/build the graph
    //--------------------------------------------------------------------------
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;

    // using -r will build a random graph
    bool random = (strcmp (matrix_name, "-r") == 0) ;

    if (!random) {
        LG_TRY (readproblem (&G, NULL,
            true, true, false, GrB_FP64, false, argc, argv)) ;
    } else {
        GrB_Index n = (argc > 2 ? atoi (argv [2]) : DEFAULT_SIZE) ;
        double density = (argc > 3 ? atof (argv [3]) : DEFAULT_DENSITY) ;
        uint64_t seed = (argc > 4 ? atoll (argv [4]) : DEFAULT_SEED) ;

        LG_TRY (LAGraph_Random_Matrix (&A, GrB_FP64, n, n, density, seed, msg)) ;
        GRB_TRY (GrB_eWiseAdd (A, NULL, NULL, GrB_PLUS_FP64, A, A, GrB_DESC_T1)) ;
        
        LG_TRY (LAGraph_New (&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg)) ;
        LG_TRY (LAGraph_Cached_NSelfEdges (G, msg)) ;
        LG_TRY (LAGraph_DeleteSelfEdges (G, msg)) ;
        // LG_TRY (LAGraph_Graph_Print (G, LAGraph_COMPLETE, stdout, msg)) ;
    }
    if (burble) {
        printf("================ DONE WITH MATRIX BUILDING ================\n") ;
    }
    GrB_Index n ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    int nt = NTHREAD_LIST ;
    
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max, nthreads_outer, nthreads_inner ;
    LG_TRY (LAGraph_GetNumThreads (&nthreads_outer, &nthreads_inner, NULL)) ;
#ifdef VERBOSE
    printf("nthreads_outer: %d, nthreads_inner: %d\n", nthreads_outer, nthreads_inner);
#endif
    nthreads_max = nthreads_outer * nthreads_inner ;
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
#ifdef VERBOSE
    printf ("threads to test: ") ;
#endif
    for (int t = 1 ; t <= nt ; t++)
    {
        int nthreads = Nthreads [t] ;
        if (nthreads > nthreads_max) continue ;
#ifdef VERBOSE
        printf (" %d", nthreads) ;
#endif
    }
#ifdef VERBOSE
    printf ("\n") ;
#endif
    if (burble) {
        printf("================ STARTING WARMUP ================\n") ;
    }
    // warmup for more accurate timing
    double tt = LAGraph_WallClockTime ( ) ;
    // GRB_TRY (LAGraph_Matrix_Print (E, LAGraph_COMPLETE, stdout, msg)) ;
    LG_TRY (LAGraph_Coarsen_Matching (&coarsened, &parent_result, &newlabel_result, &inv_newlabel_result, G, LAGraph_Matching_heavy, 0, 1, DEFAULT_SEED, msg)) ;

    tt = LAGraph_WallClockTime ( ) - tt ;

#ifdef SHOW_RESULTS
    printf("printing coarsened adjacency:\n") ;
    LG_TRY (LAGraph_Matrix_Print (coarsened, LAGraph_COMPLETE, stdout, msg)) ;
    printf("printing parent vec:\n") ;
    LG_TRY (LAGraph_Vector_Print (parent_result, LAGraph_COMPLETE, stdout, msg)) ;
    printf("printing newlabel vec:\n") ;
    LG_TRY (LAGraph_Vector_Print (newlabel_result, LAGraph_COMPLETE, stdout, msg)) ;
#endif

    GRB_TRY (GrB_free (&coarsened)) ;
    GRB_TRY (GrB_free (&parent_result)) ;
    GRB_TRY (GrB_free (&newlabel_result)) ;
    GRB_TRY (GrB_free (&inv_newlabel_result)) ;

    if (burble) {
        printf("================ WARMUP DONE ================\n") ;
    }
#ifdef VERBOSE
    printf ("warmup time %g sec\n", tt) ;
#endif

    // the GAP benchmark requires 16 trials
    int ntrials = 16 ;
    // ntrials = 1 ;    // HACK to run just one trial
#ifdef VERBOSE
    printf ("# of trials: %d\n", ntrials) ;
#endif

    for (int kk = 1 ; kk <= nt ; kk++)
    {
        int nthreads = Nthreads [kk] ;
        if (nthreads > nthreads_max) continue ;
        LG_TRY (LAGraph_SetNumThreads (1, nthreads, msg)) ;

#ifdef VERBOSE
        printf ("\n--------------------------- nthreads: %2d\n", nthreads) ;
#endif

        double total_time = 0 ;

        for (int trial = 0 ; trial < ntrials ; trial++)
        {
            int64_t seed = trial * n + 1 ;
            double tt = LAGraph_WallClockTime ( ) ;

            LG_TRY (LAGraph_Coarsen_Matching (&coarsened, &parent_result, &newlabel_result, &inv_newlabel_result, G, LAGraph_Matching_heavy, 0, 1, DEFAULT_SEED, msg)) ;

            tt = LAGraph_WallClockTime ( ) - tt ;

            GRB_TRY (GrB_free (&coarsened)) ;
            GRB_TRY (GrB_free (&parent_result)) ;
            GRB_TRY (GrB_free (&newlabel_result)) ;
            GRB_TRY (GrB_free (&inv_newlabel_result)) ;
            
#ifdef VERBOSE
            printf ("trial: %2d time: %10.7f sec\n", trial, tt) ;
#endif
            total_time += tt ;
        }

        double t = total_time / ntrials ;

#ifndef VERBOSE
        printf("%.7f\n", t);
#endif

#ifdef VERBOSE
        printf ("single-level coarsening (heavy, nopreserve, combine): %3d threads: avg time: %10.7f (sec) matrix: %s\n",
                nthreads, t, (random ? "random" : matrix_name)) ;
        fprintf (stderr, "single-level coarsening (heavy, nopreserve, combine): %3d threads: avg time: %10.7f (sec) matrix: %s\n",
                nthreads, t, (random ? "random" : matrix_name)) ;
#endif
    }

    LG_FREE_ALL ;

    LG_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
