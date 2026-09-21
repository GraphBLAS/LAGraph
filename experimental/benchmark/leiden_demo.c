//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/leiden_demo.c:
// run Leiden on a sanitized input graph and report modularity
//------------------------------------------------------------------------------

// LAGraph, (c) 2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

// Usage:
//   leiden_demo < matrixmarketfile.mtx
//   leiden_demo matrixmarketfile.mtx
//   leiden_demo matrixmarketfile.grb

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

#define LG_FREE_ALL             \
{                               \
    LAGraph_Delete (&G, NULL) ; \
    GrB_free (&c) ;             \
}

int main (int argc, char **argv)
{
#if LAGRAPH_SUITESPARSE

    char msg [LAGRAPH_MSG_LEN] ;
    msg [0] = '\0' ;
    LAGraph_Graph G = NULL ;
    GrB_Vector c = NULL ;

    bool burble = false ;
    LAGRAPH_TRY (demo_init (burble)) ;

    // read and sanitize:
    //  - make symmetric (A = A + A')
    //  - remove self edges
    //  - keep numeric values (not structural)
    //  - cast to FP64
    //  - ensure all values are nonnegative
    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, false, GrB_FP64, true, argc, argv)) ;

    LAGRAPH_TRY (LAGraph_Cached_IsSymmetricStructure (G, msg)) ;
    LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;
    LAGRAPH_TRY (LAGraph_Cached_EMin (G, msg)) ;

    double emin = 0 ;
    GRB_TRY (GrB_Scalar_extractElement_FP64 (&emin, G->emin)) ;
    if (emin < 0)
    {
        printf ("error: sanitized graph has negative edge weights (emin=%g)\n",
            emin) ;
        LG_FREE_ALL ;
        LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
        return (GrB_INVALID_VALUE) ;
    }

    uint64_t seed = 0 ;
    double t = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGraph_Leiden (&c, G, seed, msg)) ;
    t = LAGraph_WallClockTime ( ) - t ;

    double Q = 0 ;
    LAGRAPH_TRY (LAGr_Modularity (&Q, 1.0, c, G, msg)) ;

    GrB_Index n = 0, nvals = 0 ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Vector_nvals (&nvals, c)) ;

    int64_t max_label = -1 ;
    GRB_TRY (GrB_Vector_reduce_INT64 (
        &max_label, NULL, GxB_MAX_INT64_MONOID, c, NULL)) ;
    GrB_Index n_communities = (max_label < 0) ? 0 : (GrB_Index) (max_label + 1) ;

    printf ("Leiden complete\n") ;
    printf ("  n                 : %llu\n", (unsigned long long) n) ;
    printf ("  labeled nodes      : %llu\n", (unsigned long long) nvals) ;
    printf ("  communities        : %llu\n", (unsigned long long) n_communities) ;
    printf ("  modularity (gamma=1): %.12g\n", Q) ;
    printf ("  run time (sec)     : %.6g\n", t) ;

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}

