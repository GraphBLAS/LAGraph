//------------------------------------------------------------------------------
// LG_check_rcc: A hand coded RCC algorithm for CSR matircies.
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel Gomez, Texas A&M University

//------------------------------------------------------------------------------

#include <stdlib.h>
#include "LG_internal.h"
#include "LG_test.h"
#include "LG_Xtest.h"

#undef LG_FREE_WORK
#undef LG_FREE_ALL

#define LG_FREE_WORK                                    \
{                                                       \
    /* free any workspace used here */                  \
    LAGraph_Free(&a_space, NULL) ;                  \
    GrB_free (&cont) ;                                  \
    LAGraph_Free((void **)&Ai, NULL) ;               \
    LAGraph_Free((void **)&Ap, NULL) ;               \
    LAGraph_Free((void **)&slice, NULL) ;               \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    LAGraph_Free((void **)&rcc, NULL) ;     \
    GrB_free (rccs) ;      \
}
#define TIMINGS
#ifdef TIMINGS
static void print_timings (const double timings [16])
{
    double total = timings [0] + timings [1] + timings [2] + timings [3] + timings [4];
    printf ("RCC %12.6f (%4.1f%%) init\n", timings [0], 100. * timings [0] / total) ;
    printf ("RCC %12.6f (%4.1f%%) counting edges\n", timings [1], 100. * timings [1] / total) ;
    printf ("RCC %12.6f (%4.1f%%) counting nodes\n", timings [2], 100. * timings [2] / total) ;
    printf ("RCC %12.6f (%4.1f%%) cumulative sum\n", timings [3], 100. * timings [3] / total) ;
    printf ("RCC %12.6f (%4.1f%%) calculation\n", timings [4], 100. * timings [4] / total) ;
}
#endif

//Scuffed upperbound function 
static int64_t LG_binary_search    // returns upperbound - 1
(
    const int64_t pivot,
    const int64_t *LG_RESTRICT X_0,         // search in X [p_start..p_end_-1]
    const int64_t p_start,
    const int64_t p_end
)
{

    //--------------------------------------------------------------------------
    // find where the Pivot appears in X
    //--------------------------------------------------------------------------

    // binary search of X [p_start...p_end-1] for the Pivot
    int64_t pleft = p_start ;
    int64_t pright = p_end;
    while (pleft < pright)
    {
        int64_t pmiddle = pleft + (pright - pleft) / 2 ;
        bool less = (X_0 [pmiddle] < pivot) ;
        pleft  = less ? pmiddle + 1 : pleft ;
        pright = less ? pright : pmiddle ;
    }
    if(X_0[pleft] <= pivot)
        pleft++;
    return (--pleft) ;
}


int LG_check_rcc

(
    // output:
    //rccs(i): rich club coefficent of i
    GrB_Vector *rccs,    

    // input: 
    LAGraph_Graph G, //input graph
    char *msg
)
{
    #if LG_SUITESPARSE_GRAPHBLAS_V10
    GxB_Container cont = NULL;
    GrB_Matrix A = G->A;
    int64_t  *Ap = NULL, *Ai = NULL;
    void *a_space = NULL;
    GrB_Type p_type = NULL, i_type = NULL;
    int p_hand = 0, i_hand = 0;
    int n_threads = LG_nthreads_outer * LG_nthreads_inner;
    uint64_t p_n = 0, i_n = 0, p_size = 0, i_size = 0 ;
    int64_t max_deg = 0;
    uint64_t *epd = NULL, *vpd = NULL;
    int64_t *LG_RESTRICT slice  = NULL;
    double *rcc = NULL;
    #ifdef TIMINGS
    double timings [16] ;
    memset(timings, 0, 16*sizeof(double)) ;
    double tic = LAGraph_WallClockTime ( ) ;
    LG_SET_BURBLE (false) ;
    #endif

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rccs != NULL, GrB_NULL_POINTER);

    LG_ASSERT_MSG(
        G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE, 
        "G->A must be symmetric") ;
    LG_ASSERT_MSG(
        G->is_symmetric_structure == LAGraph_TRUE, GrB_INVALID_VALUE, 
        "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,
        "G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE, 
        "G->nself_edges must be zero") ; 
    GRB_TRY(GxB_Container_new(&cont)) ;
    GRB_TRY(GxB_unload_Matrix_into_Container(A, cont, NULL)) ;
    LG_ASSERT_MSG(cont->format == GxB_SPARSE, GrB_NOT_IMPLEMENTED, 
        "Matrix must be sparse") ;    
    LG_TRY (LAGraph_Malloc(
        (void **) &Ap, cont->nrows + 1, sizeof(uint64_t), NULL)) ;
    LG_TRY (LAGraph_Malloc(
        (void **) &Ai, cont->nvals, sizeof(uint64_t), NULL)) ;
    p_n = cont->nrows + 1; i_n = cont->nvals;
    GRB_TRY (GrB_Vector_extractTuples_INT64(
        NULL, Ap, &p_n, cont->p)) ;
    GRB_TRY (GrB_Vector_extractTuples_INT64(
        NULL, Ai, &i_n, cont->i)) ;
    GRB_TRY (GxB_load_Matrix_from_Container(A, cont, NULL)) ;
    GRB_TRY (GrB_Vector_reduce_INT64(
        &max_deg, NULL, GrB_MAX_MONOID_INT64, G->out_degree, NULL)) ;
    int64_t i = 0;
    #ifdef TIMINGS
    timings[0] = LAGraph_WallClockTime ( );
    #endif
    LG_TRY (LAGraph_Calloc(&a_space, max_deg * 2, sizeof(uint64_t), NULL)) ;
    LG_TRY (LAGraph_Malloc((void **)&slice, n_threads + 1, sizeof(int64_t), NULL)) ;
    epd = a_space ;
    vpd = ((uint64_t *) a_space) + max_deg ;
    LG_TRY (LAGraph_Malloc((void **) &rcc, max_deg, sizeof(double), NULL)) ;
    LG_eslice (slice, i_n, n_threads) ;
    int tid;
    #pragma omp parallel for num_threads(n_threads) schedule(static, 1) private(i)
    for (tid = 0 ; tid < n_threads ; tid++)
    {
        int64_t loc_sum = 0, dp = 0;
        int64_t loc_arr[1024];
        memset(loc_arr, 0, 1024 * sizeof(int64_t));
        i = slice[tid];
        int64_t ptr = LG_binary_search(i, Ap, 0, p_n - 1) ;
        while(i < slice[tid + 1])
        {
            while(Ap[ptr + 1] <= i) ++ptr;
            int64_t dp = Ap[ptr + 1] - Ap[ptr];
            if(dp <= 1024)
                for(; i < slice[tid + 1] && i < Ap[ptr + 1]; ++i)
                {
                    uint64_t di = Ap[Ai[i] + 1] - Ap[Ai[i]];
                    loc_arr[dp - 1] += (dp < di) + (dp <= di);
                }
            else
            {
                loc_sum = 0;
                for(; i < slice[tid + 1] && i < Ap[ptr + 1]; ++i)
                {
                    uint64_t di = Ap[Ai[i]+1] - Ap[Ai[i]];
                    loc_sum += (dp < di) + (dp <= di);
                }
                #pragma omp atomic
                    epd[dp - 1] += loc_sum ;
            }
        }
        #pragma omp critical
        {
            for(int64_t j = 0; j < 1024 && j < max_deg; ++j)
            {
                epd[j] += loc_arr[j];
            }
        }
    }
    #ifdef TIMINGS
    timings[1] = LAGraph_WallClockTime ( );
    #endif
    
    #pragma omp parallel
    {
        int64_t loc_arr[1024];
        memset(loc_arr, 0, 1024 * sizeof(int64_t));
        #pragma omp for schedule(static)
        for(i = 0; i < p_n - 1; ++i)
        {
            int64_t dp = Ap[i + 1] - Ap[i] - 1;
            if(dp < 0) continue;
            if(dp < 1024)
            {
                ++loc_arr[dp];
            }
            else
            {
                #pragma omp atomic
                    ++vpd[dp];
            }
        }  
        #pragma omp critical
        {
            for(int64_t j = 0; j < 1024 && j < max_deg; ++j)
            {
                vpd[j] += loc_arr[j];
            }
        }
    }
    
    #ifdef TIMINGS
    timings[2] = LAGraph_WallClockTime ( );
    #endif
    //run a cummulative sum (backwards)
    for(i = max_deg - 1; i > 0; --i)
    {
        vpd[i-1] += vpd[i] ;
        epd[i-1] += epd[i] ;
    }
    #ifdef TIMINGS
    timings[3] = LAGraph_WallClockTime ( );
    #endif
    #pragma omp parallel for schedule(static)
    for(i = 0; i < max_deg; ++i)
    {
        rcc[i] = ((double)epd[i]) / ((double)vpd[i] * ((double) vpd[i] - 1.0)) ;
    }
    #ifdef TIMINGS
    timings[4] = LAGraph_WallClockTime ( );
    timings[4] -= timings[3];
    timings[3] -= timings[2];
    timings[2] -= timings[1];
    timings[1] -= timings[0];
    timings[0] -= tic;
    
    print_timings(timings);
    LG_SET_BURBLE(false);
    #endif
    epd = vpd = NULL;
    GRB_TRY (GrB_Vector_new(rccs, GrB_FP64, max_deg));
    GRB_TRY (GxB_Vector_load(
        *rccs, (void **) &rcc, GrB_FP64, max_deg, max_deg * sizeof(double), 
        GrB_DEFAULT, NULL)) ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;  
    #else
    printf ("RCC not implemented for GraphBLAS versions under 10\n") ;  
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}
#undef TIMINGS
