//------------------------------------------------------------------------------
// LAGraph_BF_basic_pushpull: Bellman-Ford method for single source shortest
// paths
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_BF_basic_pushpull: Bellman-Ford single source shortest paths,
// returning just the shortest path lengths.  Contributed by Jinhao Chen and
// Tim Davis, Texas A&M.

// LAGraph_BF_basic_pushpull performs a Bellman-Ford to find out shortest path
// length from given source vertex s in the range of [0, n) on graph given as
// matrix A with size n by n. The sparse matrix A has entry A(i, j) if there is
// edge from vertex i to vertex j with weight w, then A(i, j) = w. Furthermore,
// LAGraph_BF_basic requires A(i, i) = 0 for all 0 <= i < n.

// TODO: think about the retrun values
// LAGraph_BF_basic returns GrB_SUCCESS regardless of existence of
// negative-weight cycle. However, the GrB_Vector d(k) (i.e., *pd_output) will
// be NULL when negative-weight cycle detected. Otherwise, the vector d has
// d(k) as the shortest distance from s to k.

//------------------------------------------------------------------------------

//#include "LAGraph_internal.h"
#include <LAGraph.h>
#include <LAGraphX.h>
#include <LG_internal.h>  // from src/utility

#define LAGRAPH_FREE_ALL   \
{                          \
    GrB_free(&d) ;         \
    GrB_free(&dtmp) ;      \
}


// Given a n-by-n adjacency matrix A and a source vertex s.
// If there is no negative-weight cycle reachable from s, return the distances
// of shortest paths from s as vector d. Otherwise, return d=NULL if there is
// negative-weight cycle.
// pd_output = &d, where d is a GrB_Vector with d(k) as the shortest distance
// from s to k when no negative-weight cycle detected, otherwise, d = NULL.
// A has zeros on diagonal and weights on corresponding entries of edges
// s is given index for source vertex
GrB_Info LAGraph_BF_basic_pushpull
(
    GrB_Vector *pd_output,      //the pointer to the vector of distance
    const GrB_Matrix A,         //matrix for the graph
    const GrB_Matrix AT,        //transpose of A (optional)
    const GrB_Index s           //given index of the source
)
{
    GrB_Info info;
    GrB_Index nrows, ncols, nvalA;
    // tmp vector to store distance vector after n (i.e., V) loops
    GrB_Vector d = NULL, dtmp = NULL;

    if ((A == NULL && AT == NULL) || pd_output == NULL)
    {
        // required argument is missing
        LAGRAPH_ERROR ("required arguments are NULL", GrB_NULL_POINTER) ;
    }

    *pd_output = NULL;
    bool use_vxm_with_A;
    if (A == NULL)
    {
        LAGRAPH_OK (GrB_Matrix_nrows (&nrows, AT)) ;
        LAGRAPH_OK (GrB_Matrix_ncols (&ncols, AT)) ;
        LAGRAPH_OK (GrB_Matrix_nvals (&nvalA, AT)) ;
        use_vxm_with_A = false;
    }
    else
    {
        LAGRAPH_OK (GrB_Matrix_nrows (&nrows, A)) ;
        LAGRAPH_OK (GrB_Matrix_ncols (&ncols, A)) ;
        LAGRAPH_OK (GrB_Matrix_nvals (&nvalA, A)) ;
        use_vxm_with_A = true;
    }

    // push/pull requires both A and AT
    bool push_pull = (A != NULL && AT != NULL) ;

    if (nrows != ncols)
    {
        // A must be square
        LAGRAPH_ERROR ("A must be square", GrB_INVALID_VALUE) ;
    }
    GrB_Index n = nrows;           // n = # of vertices in graph
    // average node degree
    double dA = (n == 0) ? 0 : (((double) nvalA) / (double) n) ;

    if (s >= n || s < 0)
    {
        LAGRAPH_ERROR ("invalid value for source vertex s", GrB_INVALID_VALUE) ;
    }

    // values used to determine if d should be converted to dense
    // dthreshold is used when only A or AT is available
    int64_t dthreshold ;
    if (A == NULL)
    {
        dthreshold =  LAGraph_MAX (256, sqrt ((double) n)) ;
    }
    else
    {
        dthreshold = n/2;
    }
    // refer to GraphBLAS/Source/GB_AxB_select.c
    // convert d to dense when GB_AxB_select intend to use Gustavson
    size_t csize = sizeof(double) + sizeof(int64_t) ;
    double gs_memory = ((double) n * (double) csize)/1e9 ;

    bool dsparse = true;

    // Initialize distance vector, change the d[s] to 0
    LAGRAPH_OK (GrB_Vector_new(&d, GrB_FP64, n));
    LAGRAPH_OK (GrB_Vector_setElement_FP64(d, 0, s));
    // copy d to dtmp in order to create a same size of vector
    LAGRAPH_OK (GrB_Vector_dup(&dtmp, d));
//    GxB_Vector_fprint(A, "---- A ------", GxB_SHORT, stderr);
//    GxB_Vector_fprint(d, "---- d ------", GxB_SHORT, stderr);

    int64_t iter = 0;      //number of iterations
    bool same = false;     //variable indicating if d=dtmp

    // terminate when no new path is found or more than n-1 loops
    while (!same && iter < n - 1)
    {

        double tic [2] ;
        LAGraph_Tic (tic, NULL);

        // excute semiring on d and A, and save the result to d
        if (!use_vxm_with_A)
        {
            LAGRAPH_OK (GrB_mxv(dtmp, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP64,
                AT, d, GrB_NULL));
        }
        else
        {
            LAGRAPH_OK (GrB_vxm(dtmp, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP64,
                d, A, GrB_NULL));
        }
//        double t1; LAGraph_Toc (&t1, tic, NULL) ;
        LAGRAPH_OK (LAGraph_Vector_IsEqual_type(&same, dtmp, d, GrB_FP64, NULL));
        if (!same)
        {
            GrB_Vector ttmp = dtmp;
            dtmp = d;
            d = ttmp;
//            GxB_Vector_fprint(d, "---- d ------", GxB_SHORT, stderr);
        }
        iter++;

        double t2;
        LAGraph_Toc (&t2, tic, NULL) ;
        GrB_Index dnz ;
        LAGRAPH_OK (GrB_Vector_nvals (&dnz, d)) ;

//      printf ("step %3d time1 %16.4f sec time2 %16.4f sec, nvals %.16g\n", iter, t1, t2-t1, (double) dnz) ;
//      printf ("step %3d time %16.4f sec, nvals %.16g\n", iter, t2, (double) dnz) ;
        fflush (stdout) ;

        if (dsparse)
        {
            if (!push_pull)
            {
                if (dnz > dthreshold)
                {
                    dsparse = false;
                }
            }
            else
            {
                double heap_memory = (((double) dnz+1) *
                                     5 * (double) (sizeof(int64_t))) / 1e9;
                int log2dnz = 0 ;
                while (dnz > 0)
                {
                    dnz = dnz / 2 ;
                    log2dnz++ ;
                }
                dsparse = (4 * log2dnz * heap_memory < gs_memory);
                use_vxm_with_A = dsparse;
            }

            if (!dsparse)
            {
                LAGRAPH_OK (GrB_Vector_setElement_FP64(d, 1e-16, s));
                LAGRAPH_OK (GrB_assign (d, d, NULL, INFINITY, GrB_ALL, n,
                    GrB_DESC_C)) ;
                LAGRAPH_OK (GrB_Vector_setElement_FP64(d, 0, s));
//                GxB_Vector_fprint(d, "---- d ------", GxB_SHORT, stderr);
            }
        }
    }

    // check for negative-weight cycle only when there was a new path in the
    // last loop, otherwise, there can't be a negative-weight cycle.
    if (!same)
    {
        // excute semiring again to check for negative-weight cycle
        if (!use_vxm_with_A)
        {
            LAGRAPH_OK (GrB_mxv(dtmp, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP64,
                AT, d, GrB_NULL));
        }
        else
        {
            LAGRAPH_OK (GrB_vxm(dtmp, GrB_NULL, GrB_NULL, GrB_MIN_PLUS_SEMIRING_FP64,
                d, A, GrB_NULL));
        }
        LAGRAPH_OK (LAGraph_Vector_IsEqual_type(&same, dtmp, d, GrB_FP64, NULL));

        // if d != dtmp, then there is a negative-weight cycle in the graph
        if (!same)
        {
//            GxB_Vector_fprint(d, "---- d ------", GxB_SHORT, stderr);
            // printf("A negative-weight cycle found. \n");
            LAGRAPH_FREE_ALL;
            return (GrB_SUCCESS) ;
        }
    }

    //--------------------------------------------------------------------------
    // TODO: make d sparse
    //--------------------------------------------------------------------------
    /*if (!dsparse)
    {
        LAGRAPH_OK (GrB_assign (d, d, NULL, d, GrB_ALL, n, GrB_DESC_R)) ;
        LAGRAPH_OK (GrB_Vector_setElement_FP64(d, 0, s));
        GrB_Index dnz ;
        LAGRAPH_OK (GrB_Vector_nvals (&dnz, d)) ;
        printf ("final nvals %.16g\n", (double) dnz) ;
    }*/

    (*pd_output) = d;
    d = NULL;
    LAGRAPH_FREE_ALL;
    return (GrB_SUCCESS) ;
}
