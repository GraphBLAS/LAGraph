//------------------------------------------------------------------------------
// LAGraph_BreadthFirstSearch:  breadth-first search dispatch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Breadth-first-search via push/pull method if using SuiteSparse:GraphBLAS
// and its GxB extensions, or a push-only method otherwise.  The former is
// much faster.

#include "LG_alg_internal.h"

int LAGraph_BreadthFirstSearch
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
)
{

#if LG_SUITESPARSE
    return LG_BreadthFirstSearch_SSGrB  (level, parent, G, src, pushpull, msg);
#else
    return LG_BreadthFirstSearch_vanilla(level, parent, G, src, pushpull, msg);
#endif
}

