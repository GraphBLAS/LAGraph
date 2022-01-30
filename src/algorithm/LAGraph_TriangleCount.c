//------------------------------------------------------------------------------
// LAGraph_TriangleCount:  triangle counting dispatch, basic API
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------
#define LAGraph_FREE_ALL ;

#include <LAGraph.h>
#include "LG_internal.h"
#include "LG_alg_internal.h"

//****************************************************************************
// Pick the default method with auto presort
// Compute G->ndiag, and G->rowdegree if needed.  Determine if G->A is
// symmetric, if not known.
int LAGraph_TriangleCount
(
    uint64_t *ntriangles,
    LAGraph_Graph G,
    char *msg
)
{
    // find out if graph is symmetric
    LG_TRY ( LAGraph_Property_ASymmetricStructure(G, msg) );
    LG_TRY ( LAGraph_Property_RowDegree(G, msg) );
    LG_TRY ( LAGraph_Property_NDiag(G, msg) );

    // default method and auto selection of sort
    int method  = LAGraph_TriangleCount_Default ;
    int presort = LAGraph_TriangleCount_AutoSort ;
    return LAGraph_TriangleCount_Methods (ntriangles, G, method, &presort, msg);
}

