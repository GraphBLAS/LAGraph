//------------------------------------------------------------------------------
// LG_alg_internal.h: include file for use within LAGraph algorithms
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// These definitions are not meant for the end-user of LAGraph or GraphBLAS.
// None of these methods are user-callable.

#ifndef LG_ALG_INTERNAL_H
#define LG_ALG_INTERNAL_H

#include "LG_internal.h"

int LG_BreadthFirstSearch_SSGrB
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    const LAGraph_Graph G,
    GrB_Index      src,
    char          *msg
) ;

int LG_BreadthFirstSearch_vanilla
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    const LAGraph_Graph G,
    GrB_Index      src,
    char          *msg
) ;

int LG_CC_FastSV6           // SuiteSparse:GraphBLAS method, with GxB extensions
(
    // output:
    GrB_Vector *component,  // output: array of component identifiers
    // input:
    LAGraph_Graph G,        // input graph (modified then restored)
    char *msg
) ;

int LG_CC_Boruvka
(
    // output:
    GrB_Vector *component,  // output: array of component identifiers
    // input:
    const LAGraph_Graph G,  // input graph, not modified
    char *msg
) ;

#endif
