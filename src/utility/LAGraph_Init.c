//------------------------------------------------------------------------------
// LAGraph_Init: start GraphBLAS and LAGraph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Init (char *msg)
{

    LG_CLEAR_MSG ;

    // use ANSI C memory allocation functions
    return (LAGr_Init (malloc, calloc, realloc, free, msg)) ;
}
