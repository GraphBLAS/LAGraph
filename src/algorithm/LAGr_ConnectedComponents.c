//------------------------------------------------------------------------------
// LAGr_ConnectedComponents:  connected components of an undirected graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->structure_is_symmetric must be known).

// Connected Components via LG_CC_FastSV6 if using SuiteSparse:GraphBLAS and
// its GxB extensions, or LG_CC_Boruvka otherwise.  The former is much faster.

#include "LG_alg_internal.h"

int LAGr_ConnectedComponents
(
    // output:
    GrB_Vector *component,  // component(i)=s if node i is in the component
                            // whose representative node is s
    // input:
    const LAGraph_Graph G,  // input graph
    char *msg
)
{

    #if LAGRAPH_SUITESPARSE
    return (LG_CC_FastSV6 (component, G, msg)) ;
    #else
    return (LG_CC_Boruvka (component, G, msg)) ;
    #endif
}

