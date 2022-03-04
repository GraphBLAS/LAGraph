//------------------------------------------------------------------------------
// LAGr_BreadthFirstSearch:  breadth-first search dispatch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Scott McMillan, SEI Carnegie Mellon University

//------------------------------------------------------------------------------

// Breadth-first-search via push/pull method if using SuiteSparse:GraphBLAS
// and its GxB extensions, or a push-only method otherwise.  The former is
// much faster.

// This is an Advanced algorithm (G->AT and G->rowdgree are required).

#include "LG_alg_internal.h"

int LAGr_BreadthFirstSearch
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    const LAGraph_Graph G,
    GrB_Index      src,
    char          *msg
)
{

#if LAGRAPH_SUITESPARSE
    // requires G->AT and G->rowdegree
    return LG_BreadthFirstSearch_SSGrB   (level, parent, G, src, msg) ;
#else
    // requires no properties, but G is input-only anyway
    return LG_BreadthFirstSearch_vanilla (level, parent, G, src, msg) ;
#endif
}

