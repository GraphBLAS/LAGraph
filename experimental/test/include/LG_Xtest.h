//------------------------------------------------------------------------------
// LG_Xtest.h: include file for LAGraphX test library
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#ifndef LG_XTEST_H
#define LG_XTEST_H

#include "LAGraphX.h"

int LG_check_mis        // check if iset is a valid MIS of A
(
    GrB_Matrix A,
    GrB_Vector iset,
    GrB_Vector ignore_node,     // if NULL, no nodes are ignored.  otherwise,
                        // ignore_node(i)=true if node i is to be ignored, and
                        // not added to the independent set.
    char *msg
) ;

int LG_check_ktruss
(
    // output
    GrB_Matrix *C_handle,   // the ktruss of G->A, of type GrB_UINT32
    // input
    LAGraph_Graph G,        // the structure of G->A must be symmetric
    uint32_t k,
    char *msg
) ;

#endif
