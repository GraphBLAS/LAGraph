//------------------------------------------------------------------------------
// LAGraph_TriangleCount_Methods:  triangle counting dispatch, advanced API
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

//------------------------------------------------------------------------------
#undef  LAGRAPH_FREE_ALL
#define LAGRAPH_FREE_ALL                    \
{                                           \
}

//****************************************************************************
// Advanced API: Compute G->ndiag, G->A_pattern_is_symmetric, and G->rowdegree
// (if needed) before calling.
int LAGraph_TriangleCount_Methods
(
    uint64_t       *ntriangles,
    LAGraph_Graph   G,
    int             method,
    int            *presort,
    char           *msg
)
{
#if LG_SUITESPARSE
    return LG_TriangleCount_SSGrB(ntriangles, G, method, presort, msg);
#else
    return LG_TriangleCount_vanilla(ntriangles, G, method, presort, msg);
#endif
}
