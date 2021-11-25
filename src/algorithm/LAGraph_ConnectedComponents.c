//------------------------------------------------------------------------------
// LAGraph_ConnectedComponents:  connected components of an undirected graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Connected Components via LG_CC_FastSV6 if using SuiteSparse:GraphBLAS and
// its GxB extensions, or LG_CC_Boruvka otherwise.  The former is much faster.

#include "LG_alg_internal.h"

int LAGraph_ConnectedComponents
(
    // output
    GrB_Vector *component,  // output: array of component identifiers
    // inputs
    LAGraph_Graph G,        // input graph, not modified
    char *msg
)
{
    #if LG_SUITESPARSE
    return (LG_CC_FastSV6 (component, G, msg)) ;
    #else
    return (LG_CC_Boruvka (component, G, msg)) ;
    #endif
}

