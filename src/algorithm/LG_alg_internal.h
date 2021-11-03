//------------------------------------------------------------------------------
// LG_alg_internal.h: include file for use within LAGraph algorithms
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// These definitions are not meant for the end-user of LAGraph or GraphBLAS

#ifndef LG_ALG_INTERNAL_H
#define LG_ALG_INTERNAL_H

#include <LAGraph.h>

//***************************************************************************
int LG_BreadthFirstSearch_SSGrB
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
);

//***************************************************************************
int LG_BreadthFirstSearch_vanilla
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    LAGraph_Graph  G,
    GrB_Index      src,
    bool           pushpull,
    char          *msg
);



#endif
