//------------------------------------------------------------------------------
// LAGraph_BreadthFirstSearch:  breadth-first search dispatch
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include <LAGraph.h>
#include "LG_internal.h"
#include "LG_alg_internal.h"

//****************************************************************************
int LAGraph_BreadthFirstSearch
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
)
{
#if LG_SUITESPARSE
    return LG_BreadthFirstSearch_SSGrB(level, parent,
                                       G, src, pushpull, msg);
#else
    return LG_BreadthFirstSearch_vanilla(level, parent,
                                         G, src, pushpull, msg);
#endif
}
