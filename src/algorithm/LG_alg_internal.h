//------------------------------------------------------------------------------
// LG_alg_internal.h: include file for use within LAGraph algorithms
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

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

int LG_BreadthFirstSearch_SSGrB_Extended
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    const LAGraph_Graph G,
    GrB_Index      src,
    int64_t max_level,  // < 0: no limit; otherwise, stop at this level
    int64_t dest,       // < 0: no destination; otherwise, stop if dest
                        // node is reached
    bool many_expected, // if true, the result is expected to include a fair
                        // portion of the graph.  If false, the result (parent
                        // and level) is expected to be very sparse.
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

int LG_BreadthFirstSearch_vanilla_Extended
(
    // output:
    GrB_Vector    *level,
    GrB_Vector    *parent,
    // input:
    const LAGraph_Graph G,
    GrB_Index      src,
    int64_t max_level,  // < 0: no limit; otherwise, stop at this level
    int64_t dest,       // < 0: no destination; otherwise, stop if dest
                        // node is reached
    char          *msg
) ;

int LG_CC_FastSV6           // SuiteSparse:GraphBLAS method
(
    // output:
    GrB_Vector *component,  // output: array of component identifiers
    // input:
    LAGraph_Graph G,        // input graph (modified then restored)
    char *msg
) ;

int LG_CC_FastSV7           // SuiteSparse:GraphBLAS method, with GraphBLAS v10
(
    // output:
    GrB_Vector *component,  // component(i)=r if node is in the component r
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
