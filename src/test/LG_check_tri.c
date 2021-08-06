//------------------------------------------------------------------------------
// LG_check_tri: compute the number of triangles in a graph (simple method)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// A very slow, bare-bones triangle count using a sequential saxpy-based
// method.  Computes the sum(sum((A*A).*A)), in MATLAB notation, where A is
// symmetric and treated as binary (only the pattern is used).  Diagonal
// entries are ignored.  In GraphBLAS notation, C{A} = A*A followed by
// reduce(C) to scalar.  This method is for testing only, to check the result
// of other, faster methods.  Do not benchmark this method; it is slow and
// simple by design.

#define LAGRAPH_FREE_WORK                   \
{                                           \
    LAGraph_Free ((void **) &Mark) ;        \
}

#define LAGRAPH_FREE_ALL                    \
{                                           \
    LAGRAPH_FREE_WORK ;                     \
    LAGraph_Free (&Ap) ;                    \
    LAGraph_Free (&Aj) ;                    \
    LAGraph_Free (&Ax) ;                    \
}

#include "LG_internal.h"
#include "LG_test.h"

int LG_check_tri        // -1 if out of memory, 0 if successful
(
    // output
    uint64_t *ntri,     // # of triangles in A
    // input
    LAGraph_Graph G,    // the pattern of G->A must be symmetric
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    #if !LG_SUITESPARSE
    LG_CHECK (false, -1006, "SuiteSparse:GraphBLAS required") ;
    #endif
    bool *restrict Mark = NULL ;
    GrB_Index *Ap = NULL, *Aj = NULL, *Ai = NULL ;
    void *Ax = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, ncols ;
    LG_CHECK (ntri == NULL, -1003, "ntri is NULL") ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1002, "graph is invalid") ;
    LG_CHECK (G->ndiag != 0, -104, "G->ndiag must be zero") ;
    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_pattern_is_symmetric == LAGRAPH_TRUE))
    {
        // the pattern of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -1005, "G->A must be symmetric") ;
    }
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    LG_CHECK (n != ncols, -1001, "A must be square") ;

    //--------------------------------------------------------------------------
    // allocate worspace
    //--------------------------------------------------------------------------

    Mark = (bool *) LAGraph_Calloc (n, sizeof (bool)) ;
    LG_CHECK (Mark == NULL, -1005, "out of memory") ;

    //--------------------------------------------------------------------------
    // unpack the matrix in CSR form (SuiteSparse:GraphBLAS required)
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    bool iso, jumbled ;
    GrB_TRY (GxB_Matrix_unpack_CSR (A, &Ap, &Aj, &Ax, &Ap_size, &Aj_size,
        &Ax_size, &iso, &jumbled, NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // compute the # of triangles (each triangle counted 6 times)
    //--------------------------------------------------------------------------

    // The comments below are written as if A were in CSC format, but it's
    // symmetric, so the CSR and CSC formats are the same.

    int64_t ntriangles = 0 ;
    Ai = Aj ;       // assume A is symmetric and in CSC format instead

    for (int64_t j = 0 ; j < n ; j++)
    {
        // scatter A(:,j) into Mark
        for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Mark [Ai [p]] = 1 ;
        }
        // compute sum(C(:,j)) where C(:,j) = (A * A(:,j)) .* Mark
        for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            const int64_t k = Ai [p] ;
            // C(:,j) += (A(:,k) * A(k,j)) .* Mark
            for (int64_t pa = Ap [k] ; pa < Ap [k+1] ; pa++)
            {
                // C(i,j) += (A(i,k) * A(k,j)) .* Mark
                ntriangles += Mark [Ai [pa]] ;
            }
        }
        // clear A(:,j) from Mark
        for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Mark [Ai [p]] = 0 ;
        }
    }

    //--------------------------------------------------------------------------
    // repack the matrix in CSR form for SuiteSparse:GraphBLAS
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    GrB_TRY (GxB_Matrix_pack_CSR (A, &Ap, &Aj, &Ax, Ap_size, Aj_size,
        Ax_size, iso, jumbled, NULL)) ;
    #endif

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_WORK ;
    (*ntri) = ntriangles / 6 ;
    return (0) ;
}

