//------------------------------------------------------------------------------
// LAGraph_allktruss.c: find all k-trusses of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_allktruss: find all k-trusses of a graph via GraphBLAS.
// Contributed by Tim Davis, Texas A&M.

// Given a symmetric graph A with no-self edges, LAGraph_allktruss finds all
// k-trusses of A.

// The optional output matrices Cset [3..kmax-1] are the k-trusses of A.  Their
// edges are a subset of A.  Each edge in C = Cset [k] is part of at least k-2
// triangles in C.  The structure of C is the adjacency matrix of the k-truss
// subgraph of A.  The edge weights of C are the support of each edge.  That
// is, C(i,j)=nt if the edge (i,j) is part of nt triangles in C.  All edges in
// C have support of at least k-2.  The total number of triangles in C is
// sum(C)/6.  The number of edges in C is nnz(C)/2.  C = Cset [k] is returned
// as symmetric with a zero-free diagonal, if Cset is not NULL on input.  The
// k-trusses are not returned if Cset is NULL.  Cset [kmax] is NULL since the
// kmax-truss is empty.

// The arrays ntris, nedges, and nstepss hold the output statistics.
// ntris   [k] = # of triangles in the k-truss
// nedges  [k] = # of edges in the k-truss
// nstepss [k] = # of steps required to compute the k-truss

// Usage: constructs k-trusses of A, for k = 3:kmax

//      GrB_Matrix_nrows (&n, A) ;
//      GrB_Matrix *Cset = LAGraph_malloc (n, sizeof (GrB_Matrix)) ;
//      int64_t *ntris   = LAGraph_malloc (n, sizeof (int64_t)) ;
//      int64_t *nedges  = LAGraph_malloc (n, sizeof (int64_t)) ;
//      int64_t *nstepss = LAGraph_malloc (n, sizeof (int64_t)) ;
//      GrB_Info info = LAGraph_allktruss (&Cset, A, &kmax,
//          ntris, nedges, nstepss) ;

// FIXME: add experimental/benchmark/ktruss_demo.c to benchmark k-truss
// and all-k-truss

#define LAGraph_FREE_ALL                        \
{                                               \
    if (!keep_all_ktrusses)                     \
    {                                           \
        for (int64_t kk = 3 ; kk <= k ; kk++)   \
        {                                       \
            GrB_free (&(Cset [kk])) ;           \
        }                                       \
    }                                           \
    GrB_free (&C) ;                             \
}

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// C = LAGraph_allktruss (A,k): find all k-trusses a graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_allktruss      // compute all k-trusses of a graph
(
    GrB_Matrix *Cset,           // size n, output k-truss subgraphs (optional)
    GrB_Matrix A,               // n-by-n adjacency matrix, A, not modified
    // output statistics; n3 = max (n,3)
    int64_t *kmax,              // smallest k where k-truss is empty
    int64_t *ntris,             // size n3, ntris [k] is #triangles in k-truss
    int64_t *nedges,            // size n3, nedges [k] is #edges in k-truss
    int64_t *nstepss,           // size n3, nstepss [k] is #steps for k-truss
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    int64_t k = 0 ;
    GrB_Matrix C = NULL ;
    bool keep_all_ktrusses = (Cset != NULL) ;
    LG_CHECK (A == NULL || nstepss == NULL || kmax == NULL || ntris == NULL ||
        nedges == NULL, GrB_NULL_POINTER, "input(s) are NULL") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    for (k = 0 ; k < 3 ; k++)
    {
        if (keep_all_ktrusses)
        {
            Cset [k] = NULL ;
        }
        ntris   [k] = 0 ;
        nedges  [k] = 0 ;
        nstepss [k] = 0 ;
    }
    (*kmax) = 0 ;
    GrB_Index n ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;

    //--------------------------------------------------------------------------
    // C{A} = A*A'
    //--------------------------------------------------------------------------

    GrB_Index nvals_last ;
    GrB_TRY (GrB_Matrix_nvals (&nvals_last, A)) ;
    GrB_TRY (GrB_Matrix_new (&C, GrB_UINT32, n, n)) ;
    GrB_TRY (GrB_mxm (C, A, NULL, LAGraph_plus_one_uint32, A, A,
        GrB_DESC_RST1)) ;
    int64_t nsteps = 1 ;

    //--------------------------------------------------------------------------
    // find all k-trusses
    //--------------------------------------------------------------------------

    for (k = 3 ; ; k++)
    {

        //----------------------------------------------------------------------
        // find the k-truss
        //----------------------------------------------------------------------

        while (true)
        {

            //------------------------------------------------------------------
            // C = C .* (C >= support)
            //------------------------------------------------------------------

            GrB_TRY (GrB_select (C, NULL, NULL, GrB_VALUEGE_UINT32, C, k-2,
                NULL)) ;

            //------------------------------------------------------------------
            // check if k-truss has been found
            //------------------------------------------------------------------

            GrB_Index nvals ;
            GrB_TRY (GrB_Matrix_nvals (&nvals, C)) ;
            if (nvals == nvals_last)
            {
                // k-truss has been found
                int64_t nt = 0 ;
                GrB_TRY (GrB_reduce (&nt, NULL, GrB_PLUS_MONOID_INT64, C,
                    NULL)) ;
                ntris   [k] = nt / 6 ;
                nedges  [k] = nvals / 2 ;
                nstepss [k] = nsteps ;
                nsteps = 0 ;
                if (nvals == 0)
                {
                    // this is the last k-truss
                    GrB_TRY (GrB_free (&C)) ;    // free last empty k-truss
                    (*kmax) = k ;
                    if (keep_all_ktrusses)
                    {
                        Cset [k] = NULL ;
                    }
                    return (GrB_SUCCESS) ;
                }
                else if (keep_all_ktrusses)
                {
                    // save the k-truss in the list of output k-trusses
                    // TODO: if Cset [k] == Cset [k-1], then do not save it.
                    // Set it to NULL to denote that the k-truss is the
                    // same as the (k-1)-truss.  Also, advance quickly to
                    // the next k, setting k = min (C).
                    GrB_TRY (GrB_Matrix_dup (&(Cset [k]), C)) ;
                }
                // start finding the next k-truss
                break ;
            }

            // continue searching for this k-truss
            nvals_last = nvals ;
            nsteps++ ;

            //------------------------------------------------------------------
            // C{C} = C*C'
            //------------------------------------------------------------------

            GrB_TRY (GrB_mxm (C, C, NULL, LAGraph_plus_one_uint32, C, C,
                GrB_DESC_RST1)) ;
        }
    }
}

