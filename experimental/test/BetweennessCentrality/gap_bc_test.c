//------------------------------------------------------------------------------
// bc_gap_test: betweenness centrality for the GAP benchmark
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

#define NTHREAD_LIST 1
// #define NTHREAD_LIST 2
#define THREAD_LIST 0

// #define NTHREAD_LIST 4
// #define THREAD_LIST 40, 20, 10, 1

// #define NTHREAD_LIST 6
// #define THREAD_LIST 64, 32, 24, 12, 8, 4

// Contributed by Scott Kolodziej and Tim Davis, Texas A&M University

// usage:
// bc_gap_test matrixfile.mtx sourcenodes.mtx
// in is the Matrix Market file, out is the level set.

#include "bc_test.h"
#include "../../../GraphBLAS/Source/GB_Global.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    GrB_free (&A);                  \
    GrB_free (&AT);                 \
    GrB_free (&Abool);              \
    GrB_free (&v);                  \
    GrB_free (&v_brandes);          \
    GrB_free (&v_batch);            \
    GrB_free (&v_batch4);           \
    GrB_free (&v_batch5);           \
    GrB_free (&SourceNodes) ;       \
}

int main (int argc, char **argv)
{
    GrB_Info info;

    GrB_Matrix A = NULL ;
    GrB_Matrix AT = NULL ;
    GrB_Matrix Abool = NULL ;
    GrB_Vector v = NULL ;
    GrB_Vector v_brandes = NULL ;
    GrB_Vector v_batch = NULL ;
    GrB_Vector v_batch4 = NULL ;
    GrB_Vector v_batch5 = NULL ;
    GrB_Matrix SourceNodes = NULL ;
    LAGRAPH_OK (LAGraph_init ( )) ;
    LAGRAPH_OK (GxB_set (GxB_BURBLE, false)) ;

    printf ("using: %s v%d.%d.%d [%s]\n",
        GxB_IMPLEMENTATION_NAME,
        GxB_IMPLEMENTATION_MAJOR,
        GxB_IMPLEMENTATION_MINOR,
        GxB_IMPLEMENTATION_SUB,
        GxB_IMPLEMENTATION_DATE) ;

    uint64_t seed = 1;
    bool tests_pass = true;

    int nt = NTHREAD_LIST ;
    int Nthreads [20] = { 0, THREAD_LIST } ;
    int nthreads_max = LAGraph_get_nthreads ( ) ;
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

    // GxB_set (GxB_NTHREADS, 1) ;
    // GxB_set (GxB_CHUNK, 1) ;

    // Start the timer
    double tic [2];
    LAGraph_tic (tic);

    //--------------------------------------------------------------------------
    // read in a matrix from a file and convert to boolean
    //--------------------------------------------------------------------------

    int batch_size = 4 ;
    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ; 

    if (argc > 1)
    {
        // Usage:
        //      ./bc_gap_test matrixfile.mtx sourcenodes.mtx
        //      ./bc_gap_test matrixfile.grb sourcenodes.mtx

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
            LAGRAPH_OK (LAGraph_binread (&A, filename)) ;
        }
        else
        {
            printf ("Reading Matrix Market file: %s\n", filename) ;
            FILE *f = fopen (filename, "r") ;
            if (f == NULL)
            {
                printf ("Matrix file not found: [%s]\n", filename) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_mmread(&A, f));
            fclose (f) ;
        }

        // read in source nodes in Matrix Market format from the input file
        if (argc > 2)
        {
            char *sourcefile = argv [2] ;
            printf ("sources: %s\n", sourcefile) ;
            FILE *f = fopen (sourcefile, "r") ;
            if (f == NULL)
            {
                printf ("Source node file not found: [%s]\n", sourcefile) ;
                exit (1) ;
            }
            LAGRAPH_OK (LAGraph_mmread (&SourceNodes, f)) ;
            fclose (f) ;
        }
    }
    else
    {

        // Usage:  ./bc_gap_test < matrixfile.mtx
        printf ("matrix: from stdin\n") ;

        // read in the file in Matrix Market format from stdin
        LAGRAPH_OK (LAGraph_mmread(&A, stdin));
    }

    double t_read = LAGraph_toc (tic) ;
    printf ("read time: %g sec\n", t_read) ;

    LAGraph_tic (tic);

    // convert to pattern-only
    LAGRAPH_OK (LAGraph_pattern(&Abool, A, GrB_BOOL));

    GrB_free (&A);
    A = Abool;
    Abool = NULL ;

    // finish any pending computations
    GrB_Index nvals;
    GrB_Matrix_nvals (&nvals, A);

    //--------------------------------------------------------------------------
    // get the size of the problem.
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols;
    LAGRAPH_OK (GrB_Matrix_nrows(&nrows, A));
    LAGRAPH_OK (GrB_Matrix_ncols(&ncols, A));
    GrB_Index n = nrows;

    //--------------------------------------------------------------------------
    // get the source nodes
    //--------------------------------------------------------------------------

    #define NSOURCES 32

    if (SourceNodes == NULL)
    {
        LAGRAPH_OK (GrB_Matrix_new (&SourceNodes, GrB_INT64, NSOURCES, 1)) ;
        srand (1) ;
        for (int k = 0 ; k < NSOURCES ; k++)
        {
            int64_t i = 1 + (rand ( ) % n) ;    // in range 1 to n
            // SourceNodes [k] = i 
            LAGRAPH_OK (GrB_Matrix_setElement (SourceNodes, i, k, 0)) ;
        }
    }
    GrB_Matrix_nvals (&nvals, SourceNodes);

    GrB_Index nsource ;
    LAGRAPH_OK (GrB_Matrix_nrows(&nsource, SourceNodes));
    if (nsource % batch_size != 0)
    {
        printf ("SourceNode size must be multiple of batch_size (%d)\n",
            batch_size) ;
        exit (1) ;
    }

    double t_setup = LAGraph_toc (tic) ;
    // printf ("setup time: %g sec\n", t_setup) ;

    // AT = A'
    LAGraph_tic (tic);
    bool A_is_symmetric =
        (n == 134217726 ||  // HACK for kron
         n == 134217728) ;  // HACK for urand
    if (!A_is_symmetric)
    {
        LAGRAPH_OK (GrB_Matrix_new (&AT, GrB_BOOL, n, n)) ;
        LAGRAPH_OK (GrB_transpose (AT, NULL, NULL, A, NULL)) ;
        LAGRAPH_OK (LAGraph_isequal (&A_is_symmetric, A, AT, NULL)) ;
    }
    if (A_is_symmetric)
    {
        printf ("A is symmetric\n") ;
        GrB_free (&AT) ;
    }
    else
    {
        printf ("A is unsymmetric\n") ;
    }
    double t_transpose = LAGraph_toc (tic) ;
    printf ("transpose time: %g\n", t_transpose) ;

    #define K 1024
    #define M (K*K)

    /*
    int nt = 7 ;
    int Nthreads [7+1] = { 0,
        1, 2, 4, 8, 12, 20, 40 } ;
    */

    /*
    int nt = 3 ;
    int Nthreads [6+1] = { 0,
        10, 20, 40 } ;        // hypersparse
    */

    /*
    int nt = 4 ;
    int Nthreads [6+1] = { 0,
        1, 2, 4, 8 } ;              // slash
    */

    //--------------------------------------------------------------------------
    // Begin tests
    //--------------------------------------------------------------------------

    GrB_Matrix_nvals (&nvals, A);
    printf ("\n========== input graph: nodes: %"PRIu64" edges: %"PRIu64"\n",
        n, nvals) ;

    int ntrials = 0 ;
    double total_time_1 = 0 ;
    double total_time_x3 [8+1] ;
    double total_time_4  [8+1] ;
    double total_time_5  [8+1] ;
    double total_timing  [3] = { 0, 0, 0 } ;

    for (int t = 0 ; t < 9 ; t++)
    {
        total_time_x3 [t] = 0 ;
        total_time_4  [t] = 0 ;
        total_time_5  [t] = 0 ;
    }


    #if 0
    {
        char filename [256] ;
        sprintf (filename, "batch_src_%lu.mtx", n) ;
        FILE *f = fopen (filename, "w") ;
        LAGraph_mmwrite (SourceNodes, f) ;
        fclose (f) ;
    }
    #endif

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
            LAGRAPH_OK (GrB_Matrix_extractElement (&source, SourceNodes,
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

//      LAGRAPH_OK (LAGraph_bc_batch  (&v_batch, A, vertex_list, batch_size)) ;
//      LAGRAPH_OK (LAGraphX_bc_batch (&v_batch, A, vertex_list, batch_size)) ;
//      LAGRAPH_OK (LAGraphX_bc_batch2 (&v_batch, A, vertex_list, batch_size)) ;


#if 0
        printf ("---\n") ;
        // version X3
        {
            for (int t = 1 ; t <= nt ; t++)
            {
                if (Nthreads [t] > nthreads_max) continue ;
                GxB_set (GxB_NTHREADS, Nthreads [t]) ;
                double timing [3] = { 0, 0, 0 } ;
                GrB_free (&v_batch) ;
                LAGraph_tic (tic) ;
                LAGRAPH_OK (LAGraphX_bc_batch3 (&v_batch, A,
                    ((AT == NULL) ? A : AT),
                    vertex_list, batch_size, timing)) ;
                double t2 = LAGraph_toc (tic) ;
        //      total_timing [0] += timing [0] ;        // pushpull
        //      total_timing [1] += timing [1] ;        // allpush
        //      total_timing [2] += timing [2] ;        // allpull
                total_time_x3 [t] += t2 ;
                printf ("Batch X3 time %2d: %12.4f (sec), rate: %10.3f\n",
                    Nthreads [t], t2, 1e-6*((double) nvals) / t2) ;
            }
        }
        // GxB_print (v_batch, 2) ;
        printf ("---\n") ;
#endif


#if 0
        // version 4
//      GB_Global_timing_clear_all ( ) ;
        {
            for (int t = 1 ; t <= nt ; t++)
            {
                if (Nthreads [t] > nthreads_max) continue ;
                GxB_set (GxB_NTHREADS, Nthreads [t]) ;
                GrB_free (&v_batch4) ;
                LAGraph_tic (tic) ;
                LAGRAPH_OK (LAGraph_bc_batch4 (&v_batch4, A,
                    ((AT == NULL) ? A : AT),
                    vertex_list, batch_size)) ;
                double t2 = LAGraph_toc (tic) ;
                printf ("Batch v4 time %2d: %12.4f (sec)\n",
                    Nthreads [t], t2) ;
                total_time_4 [t] += t2 ;
            }
        }
#endif

//      for (int k = 0 ; k < 20 ; k++)
//      {
//          double t = GB_Global_timing_get (k) ;
//          if (t > 0) printf ("phase %2d: %12.4f msec\n", k, t*1e3) ;
//      }

        // back to default
        GxB_set (GxB_NTHREADS, nthreads_max) ;

        // GxB_print (v_batch4, 2) ;

        // version 5
//      GB_Global_timing_clear_all ( ) ;
        {
            for (int t = 1 ; t <= nt ; t++)
            {
                if (Nthreads [t] > nthreads_max) continue ;
                GxB_set (GxB_NTHREADS, Nthreads [t]) ;
                GrB_free (&v_batch5) ;
                LAGraph_tic (tic) ;
                LAGRAPH_OK (LAGraph_bc_batch5 (&v_batch5, A,
                    ((AT == NULL) ? A : AT),
                    vertex_list, batch_size)) ;
                double t2 = LAGraph_toc (tic) ;
                printf ("Batch v5 time %2d: %12.4f (sec)\n",
                    Nthreads [t], t2) ;
                // GxB_print (v_batch5, 2) ;
                fflush (stdout) ;
                total_time_5 [t] += t2 ;
            }
        }

#if 0
        LAGRAPH_OK (GrB_Vector_new(&v_batch, GrB_FP64, n));
        for (int64_t k = 0 ; k < 4 ; k++)
        {
            // get the kth source node
            GrB_Index source = vertex_list [k] ;
//          LAGRAPH_OK (LAGraph_bc (&v, A, source)) ;
            LAGRAPH_OK (LAGraph_bc2 (&v, A, source)) ;
            LAGRAPH_OK (GrB_eWiseAdd(v_batch, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, v_batch, v, GrB_NULL));
            GrB_free (&v) ;
        }
#endif

        GrB_Type type ;
        GxB_Vector_type (&type, v_batch4) ;
        double tol ;
        if (type == GrB_FP32)
        {
            // printf("LAGraph batch is FP32\n");
            tol = 0.005 ;
        }
        else
        {
            // printf("LAGraph batch is FP64\n");
            tol = 1e-10 ;
        }

        //----------------------------------------------------------------------
        // Compute betweenness centrality from four source nodes (Brandes)
        //----------------------------------------------------------------------

// comment this out to skip the test of Brandes' method
#if 0
        // Start the timer
        LAGraph_tic (tic);

        LAGRAPH_OK (GrB_Vector_new(&v_brandes, GrB_FP64, n));
        for (int64_t k = 0 ; k < batch_size ; k++)
        {
            // get the kth source node
            GrB_Index source = vertex_list [k] ;
            LAGRAPH_OK (LAGraph_bc (&v, A, source)) ;
//          LAGRAPH_OK (LAGraph_bc2 (&v, A, source)) ;
            LAGRAPH_OK (GrB_eWiseAdd(v_brandes, GrB_NULL, GrB_NULL,
                GrB_PLUS_FP64, v_brandes, v, GrB_NULL));
            GrB_free (&v) ;
        }

        // GxB_print (v_brandes, 2) ;

        // Stop the timer
        double t1 = LAGraph_toc (tic) ;
        printf ("Brandes  time: %12.4f (sec), rate: %g (1e6 edges/sec)\n",
            t1, 1e-6*((double) nvals) / t1) ;

        total_time_1 += t1 ;
#endif

        //----------------------------------------------------------------------
        // check result
        //----------------------------------------------------------------------

        if (kstart == 0 && v_brandes != NULL)
        {
            printf ("max relative error: ") ;
            fflush (stdout) ;
            // printf ("checking result ... \n") ;
            double maxerr_x3 = 0 ;
            double maxerr_4 = 0 ;
            double xmax = 0 ;
            for (int64_t i = 0; i < n; i++)
            {
                // TODO this test is slow.  Use GrB_eWiseAdd and GrB_reduce.

                double x1 = 0;
                LAGRAPH_OK (GrB_Vector_extractElement (&x1, v_brandes, i));
                xmax = fmax (xmax, fabs (x1)) ;

                double x2 = 0;
                if (v_batch != NULL)
                {
                    LAGRAPH_OK (GrB_Vector_extractElement (&x2, v_batch, i));
                }
                double err = fabs (x1 - x2) ;
                maxerr_x3 = fmax (maxerr_x3, err) ;

                double x4 = 0;
                if (v_batch4 != NULL)
                {
                    LAGRAPH_OK (GrB_Vector_extractElement (&x4, v_batch4, i));
                }
                err = fabs (x1 - x4) ;
                maxerr_4 = fmax (maxerr_4, err) ;

                // Check that both methods give the same results
                /*
                bool test_result = err < tol ;
                tests_pass &= test_result;
                if (!test_result)
                {
                    printf ("%"PRId64 " %g %g %g ", i, x1, x2, err) ;
                    printf ("FAIL") ;
                    printf ("   Failure at index %"PRId64"\n", i);
                    printf ("   x1 = %g\n", x1);
                    printf ("   x2 = %g\n", x2);
                    printf ("   Error = %f\n", err) ;
                    printf ("\n") ;
                    break ;
                }
                */
            }
            printf ("%g %g\n", maxerr_x3 / xmax, maxerr_4 / xmax) ;

            #if 0
            printf ("writing results to mtx files:\n") ;
            FILE *f = fopen ("brandes_result.mtx", "w") ;
            LAGraph_mmwrite ((GrB_Matrix) v_brandes, f) ;
            fclose (f) ;
            #endif
        }

        #if 0
        {
            char filename [256] ;
            sprintf (filename, "batch_%02ld_%lu.mtx", kstart, n) ;
            FILE *f = fopen (filename, "w") ;
            LAGraph_mmwrite ((GrB_Matrix) v_batch4, f) ;
            fclose (f) ;
        }
        #endif

        GrB_free (&v_brandes) ;
        GrB_free (&v_batch) ;
        GrB_free (&v_batch4) ;
        GrB_free (&v_batch5) ;

        // HACK: uncomment this to just do the first batch
        // break ;
    }

    //--------------------------------------------------------------------------
    // free all workspace and finish
    //--------------------------------------------------------------------------

    printf ("\nntrials: %d\n", ntrials) ;
    if (total_time_1 > 0 && ntrials > 1)
    {
        printf ("Average time per trial (Brandes): %g sec\n",
            total_time_1 / ntrials);
    }

    if (total_time_x3 [1] > 0) // && ntrials > 1)
    {
        for (int t = 1 ; t <= nt ; t++)
        {
            if (Nthreads [t] > nthreads_max) continue ;
            double t2 = total_time_x3 [t] / ntrials ;
            printf ("Ave (BatchX3) %2d: %10.3f sec, rate %10.3f\n",
                Nthreads [t], t2, 1e-6*((double) nvals) / t2) ;
            if (n > 2000)
            {
                LAGr_log (matrix_name, "BatchX3", Nthreads [t], t2) ;
            }
        }
    }

    printf ("\n") ;

    if (total_time_4 [1] > 0) // && ntrials > 1)
    {
        printf ("\n") ;
        for (int t = 1 ; t <= nt ; t++)
        {
            if (Nthreads [t] > nthreads_max) continue ;
            double t2 = total_time_4 [t] / ntrials ;
            printf ("Ave (Batch4)  %2d: %10.3f sec, rate %10.3f\n",
                Nthreads [t], t2, 1e-6*((double) nvals) / t2) ;
            if (n > 2000)
            {
                LAGr_log (matrix_name, "Batch4", Nthreads [t], t2) ;
            }
        }
    }

    if (total_time_5 [1] > 0) // && ntrials > 1)
    {
        printf ("\n") ;
        for (int t = 1 ; t <= nt ; t++)
        {
            if (Nthreads [t] > nthreads_max) continue ;
            double t2 = total_time_5 [t] / ntrials ;
            printf ("Ave (Batch5)  %2d: %10.3f sec, rate %10.3f\n",
                Nthreads [t], t2, 1e-6*((double) nvals) / t2) ;
            fprintf (stderr, "Avg: BC (batch5)  %3d: %10.3f sec: %s\n",
                Nthreads [t], t2, matrix_name) ;
            if (n > 2000)
            {
                LAGr_log (matrix_name, "Batch5", Nthreads [t], t2) ;
            }
        }
    }

    LAGRAPH_FREE_ALL;
    LAGRAPH_OK (LAGraph_finalize());
    if (tests_pass)
    {
        printf("%s: all tests passed\n", argv [0]);
    }
    else
    {
        printf("%s: TEST FAILURE\n", argv [0]);
    }
    return (GrB_SUCCESS);
}

