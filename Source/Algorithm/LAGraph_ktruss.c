//------------------------------------------------------------------------------
// LAGraph_ktruss.c: find the k-truss subgraph of a graph via GraphBLAS
//------------------------------------------------------------------------------

// Given a symmetric graph A with no-self edges, ktruss_graphblas finds the
// k-truss subgraph of A.

// The edge weights of A are treated as binary.  Explicit zero entries in A are
// treated as non-edges.  Any type will work, but uint32 is recommended for
// fastest results since that is the type used here for the semiring.
// GraphBLAS will do typecasting internally, but that takes extra time. 

// The output matrix C is the k-truss subgraph of A.  Its edges are a subset of
// A.  Each edge in C is part of at least k-2 triangles in C.  The pattern of C
// is the adjacency matrix of the k-truss subgraph of A.  The edge weights of C
// are the support of each edge.  That is, C(i,j)=nt if the edge (i,j) is part
// of nt triangles in C.  All edges in C have support of at least k-2.  The
// total number of triangles in C is sum(C)/6.  The number of edges in C is
// nnz(C)/2.  C is returned as symmetric with a zero-free diagonal.

// Usage: constructs C as the k-truss of A
//      GrB_Matrix C = NULL ;
//      int32_t nsteps ;
//      GrB_Info info = LAGraph_ktruss (&C, A, k, &nsteps) ;

// Compare this function with the MATLAB equivalent, ktruss.m, in
// LAGraph/Test/KTruss.  This function is derived from SuiteSparse/GraphBLAS/
// Extras/ktruss_graphblas.c.

// TODO add cite

#define LAGRAPH_FREE_ALL                \
    GrB_free (&C) ;

#include "LAGraph_internal.h"

//------------------------------------------------------------------------------
// C = ktruss_graphblas (A,k): find the k-truss subgraph of a graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_ktruss         // compute the k-truss of a graph
(
    GrB_Matrix *Chandle,        // output k-truss subgraph, C
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

    if (Chandle == NULL) return (GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Index n ;
    GrB_Matrix C = NULL ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    LAGRAPH_OK (GrB_Matrix_new (&C, GrB_UINT32, n, n)) ;

    // for the select operator
    uint32_t support = (k-2) ;

    // last_cnz = nnz (A)
    GrB_Index cnz, last_cnz ;
    LAGRAPH_OK (GrB_Matrix_nvals (&last_cnz, A)) ;

    //--------------------------------------------------------------------------
    // find the k-truss of A
    //--------------------------------------------------------------------------

    for (int32_t nsteps = 1 ; ; nsteps++)
    {

        //----------------------------------------------------------------------
        // C<C> = C*C
        //----------------------------------------------------------------------

        GrB_Matrix C2 = (nsteps == 1) ? A : C ;
        LAGRAPH_OK (GrB_mxm (C, C2, NULL, GxB_PLUS_LAND_UINT32, C2, C2, NULL)) ;

        //----------------------------------------------------------------------
        // C = C .* (C >= support)
        //----------------------------------------------------------------------

        LAGRAPH_OK (GxB_select (C, NULL, NULL, LAGraph_support, C, &support,
            NULL)) ;

        //----------------------------------------------------------------------
        // check if the k-truss has been found
        //----------------------------------------------------------------------

        LAGRAPH_OK (GrB_Matrix_nvals (&cnz, C)) ;
        if (cnz == last_cnz)
        {
            (*Chandle) = C ;                        // return the result C
            if (p_nsteps != NULL)
            {
                (*p_nsteps) = nsteps ;              // return # of steps
            }
            return (GrB_SUCCESS) ;
        }
        last_cnz = cnz ;
    }
}

