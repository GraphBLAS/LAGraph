//------------------------------------------------------------------------------
// LAGraph_TriangleCount:  triangle counting dispatch, basic API
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Scott McMillan, SEI Carnegie Mellon University

//------------------------------------------------------------------------------

// This is a Basic algorithm (G->ndiag, G->row_degree, G->structure_is_symmetric
// are computed, if not present).

#define LG_FREE_ALL ;

#include <LAGraph.h>
#include "LG_internal.h"
#include "LG_alg_internal.h"

//****************************************************************************
// Pick the default method with auto presort
// Compute G->ndiag, and G->row_degree if needed.  Determine if G->A is
// symmetric, if not known.

int LAGraph_TriangleCount
(
    // output:
    uint64_t      *ntriangles,   // # of triangles
    // input/output:
    LAGraph_Graph  G,       // G->ndiag, G->row_degree, G->structure_is_symmetric
                            // are computed, if not already present
    char          *msg
)
{
    // find out if graph is symmetric, compute G->row_degree, and G->ndiag
    LG_TRY (LAGraph_Cached_SymmetricStructure (G, msg)) ;
    LG_TRY (LAGraph_Cached_RowDegree (G, msg)) ;
    LG_TRY (LAGraph_Cached_NDiag (G, msg)) ;

    // default method and auto selection of sort
    int method  = LAGraph_TriangleCount_Default ;
    int presort = LAGraph_TriangleCount_AutoSort ;
    return (LAGr_TriangleCount (ntriangles, G, method, &presort, msg)) ;
}

