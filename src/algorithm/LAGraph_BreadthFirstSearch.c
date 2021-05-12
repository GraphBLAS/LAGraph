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

//***************************************************************************
int LAGraph_BreadthFirstSearch_SSGrB
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
);

//***************************************************************************
int LAGraph_BreadthFirstSearch_vanilla
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
);

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
#if !defined(LG_VANILLA) && defined(GxB_SUITESPARSE_GRAPHBLAS)
    return LAGraph_BreadthFirstSearch_SSGrB(level, parent,
                                            G, src, pushpull, msg);
#else
    return LAGraph_BreadthFirstSearch_vanilla(level, parent,
                                              G, src, pushpull, msg);
#endif
}
