//------------------------------------------------------------------------------
// LAGraph_prune_diag: remove diagonal entries from a matrix
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
// FIXME: this is not yet included in the test coverage suite
// FIXME: remove this and use LAGraph_DeleteDiag.

// LAGraph_prune_diag: contributed by Tim Davis.  Removes all diagonal entries
// from a matrix.

//------------------------------------------------------------------------------

#define LAGraph_FREE_ALL GrB_free (&M)

#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
GrB_Info LAGraph_prune_diag // remove all entries from the diagonal
(
    GrB_Matrix A
)
{

    GrB_Info info ;
    GrB_Matrix M = NULL ;

    GrB_Index m, n ;
    LAGRAPH_OK (GrB_Matrix_nrows (&m, A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    int64_t k = LAGraph_MIN (m, n) ;

    // M = diagonal mask matrix
    LAGRAPH_OK (GrB_Matrix_new (&M, GrB_BOOL, m, n)) ;
    for (int64_t i = 0 ; i < k ; i++)
    {
        // M(i,i) = true ;
        LAGRAPH_OK (GrB_Matrix_setElement (M, (bool) true, i, i)) ;
    }

    // remove self edges (via M)
    LAGRAPH_OK (GrB_assign (A, M, NULL, A, GrB_ALL, m, GrB_ALL, n,
                            GrB_DESC_RC)) ;

    LAGraph_FREE_ALL ;
    return (GrB_SUCCESS) ;
}
