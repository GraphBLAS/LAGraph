//------------------------------------------------------------------------------
// LAGraph_ktruss: k-truss subgraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_ktruss: k-truss subgraph, contributed by Tim Davis, Texas A&M.

// Given a symmetric graph A with no-self edges, ktruss_graphblas finds the
// k-truss subgraph of A.

// TODO add sanitize step to remove diagonal

// FIXME:  currently relies on a GxB* function in SuiteSparse/GraphBLAS v3.0.1,
// including the new GxB_Scalar.  Write a version using pure GraphBLAS,
// using GrB_select

// The edge weights of A are treated as binary.  Explicit zero entries in A are
// treated as non-edges.  Any type will work, but uint32 is recommended for
// fastest results since that is the type used here for the semiring.
// GraphBLAS will do typecasting internally, but that takes extra time.

// The output matrix C is the k-truss subgraph of A.  Its edges are a subset of
// A.  Each edge in C is part of at least k-2 triangles in C.  The structure of
// C is the adjacency matrix of the k-truss subgraph of A.  The edge weights of
// C are the support of each edge.  That is, C(i,j)=nt if the edge (i,j) is
// part of nt triangles in C.  All edges in C have support of at least k-2.
// The total number of triangles in C is sum(C)/6.  The number of edges in C is
// nnz(C)/2.  C is returned as symmetric with a zero-free diagonal.

// Usage: constructs C as the k-truss of A
//      GrB_Matrix C = NULL ;
//      int32_t nsteps ;
//      GrB_Info info = LAGraph_ktruss (&C, A, k, &nsteps) ;

// Compare this function with the MATLAB equivalent, ktruss.m, in
// LAGraph/Test/KTruss.  This function is derived from SuiteSparse/GraphBLAS/
// Extras/ktruss_graphblas.c.

// TODO add cite

#define LAGraph_FREE_ALL                \
    GrB_free (&Support) ;               \
    GrB_free (&LAGraph_support) ;       \
    GrB_free (&C) ;

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************

static bool LAGraph_support_function 
(
const GrB_Index i, const GrB_Index j,
const uint32_t *x, const uint32_t *support)
{
    return ((*x) >= (*support)) ;
}

//------------------------------------------------------------------------------
// C = ktruss_graphblas (A,k): find the k-truss subgraph of a graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_ktruss         // compute the k-truss of a graph
(
    GrB_Matrix *Chandle,        // output k-truss subgraph, C
    GrB_Type *C_type,
    const GrB_Matrix A,         // input adjacency matrix, A, not modified
    const uint32_t k,           // find the k-truss, where k >= 3
    int32_t *p_nsteps           // # of steps taken (ignored if NULL)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // ensure k is 3 or more
    if (k < 3) return (GrB_INVALID_VALUE) ;

    if ((Chandle == NULL) || (C_type == NULL))
        return (GrB_NULL_POINTER) ;
    (*Chandle) = NULL ;
    (*C_type) = NULL;

#if LG_SUITESPARSE

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Index n ;
    GrB_Matrix C = NULL ;
    GxB_Scalar Support = NULL ;
    GxB_SelectOp LAGraph_support = NULL ;

    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_UINT32, n, n)) ;
    *C_type = GrB_UINT32;

    // for the select operator
    uint32_t support = (k-2) ;
    LAGRAPH_OK (GxB_Scalar_new (&Support, GrB_UINT32)) ;
    LAGRAPH_OK (GxB_Scalar_setElement (Support, support)) ;
    GrB_Index ignore ;
    LAGRAPH_OK (GxB_Scalar_nvals (&ignore, Support)) ;

    // last_cnz = nnz (A)
    GrB_Index cnz, last_cnz ;
    LAGRAPH_OK (GrB_Matrix_nvals (&last_cnz, A)) ;

    // allocate the select function for ktruss and allktruss
    LAGRAPH_OK (GxB_SelectOp_new (&LAGraph_support,
        (GxB_select_function) LAGraph_support_function,
        GrB_UINT32, GrB_UINT32)) ;

    //--------------------------------------------------------------------------
    // find the k-truss of A
    //--------------------------------------------------------------------------

    for (int32_t nsteps = 1 ; ; nsteps++)
    {

        double tic [2] ;
        LAGraph_Tic (tic, NULL) ;

        //----------------------------------------------------------------------
        // C<C> = C*C
        //----------------------------------------------------------------------

        GrB_Matrix C2 = (nsteps == 1) ? A : C ;
        LAGRAPH_OK (GrB_mxm (C, C2, NULL, GxB_PLUS_LAND_UINT32, C2, C2, NULL)) ;

        //----------------------------------------------------------------------
        // C = C .* (C >= support)
        //----------------------------------------------------------------------

        LAGRAPH_OK (GxB_select (C, NULL, NULL, LAGraph_support, C, Support,
                                NULL)) ;

        //----------------------------------------------------------------------
        // check if the k-truss has been found
        //----------------------------------------------------------------------

        LAGRAPH_OK (GrB_Matrix_nvals (&cnz, C)) ;
        double t;
        LAGraph_Toc (&t, tic, NULL) ;
        printf ("step %3d time %16.4f sec, nvals %.16g\n",
                nsteps, t, (double) cnz) ;
        fflush (stdout) ;
        if (cnz == last_cnz)
        {
            (*Chandle) = C ;                        // return the result C
            if (p_nsteps != NULL)
            {
                (*p_nsteps) = nsteps ;              // return # of steps
            }
            GrB_free (&Support) ;  Support = NULL;
            return (GrB_SUCCESS) ;
        }
        last_cnz = cnz ;
    }

#else
    // FIXME: create a pure GrB version with the v2 spec
    return (GrB_NO_VALUE) ;
#endif
}
