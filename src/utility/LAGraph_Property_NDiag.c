//------------------------------------------------------------------------------
// LAGraph_Property_NDiag: count the # of diagonal entries of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_Property_NDiag  // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to compute G->ndiag for
    char *msg
)
{
    //--------------------------------------------------------------------------
    // clear msg and check G
    //--------------------------------------------------------------------------

    LG_CHECK_INIT (G, msg) ;

    // already computed
    if (G->ndiag != LAGRAPH_UNKNOWN)
    {
        return 0;
    }

    //--------------------------------------------------------------------------
    // compute G->ndiag
    //--------------------------------------------------------------------------
    return (LG_ndiag (&G->ndiag, G->A, G->A_type, msg)) ;
}
