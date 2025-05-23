//------------------------------------------------------------------------------
// LAGraph_SwapEdges: randomly swaps edges in a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Gabriel Gomez, Texas A&M University

//------------------------------------------------------------------------------

// References:

// R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, and U. Alon, “On the 
// uniform generation of random graphs with prescribed degree sequences,” 2004.
#include "LG_internal.h"
#include "LAGraphX.h"
int LAGraph_SwapEdges
(
    // output
    LAGraph_Graph *G_new,  // A new graph with the same degree for each node
    double *pQ,            // Actual Swaps proformed per edge
    // input: not modified
    const LAGraph_Graph G, // Graph to be randomized.
    double Q,              // Swaps per edge
    char *msg
)
{
    LG_ASSERT(pQ != NULL, GrB_NULL_POINTER);
    LG_ASSERT(Q > 0.0, GrB_INVALID_VALUE);
    GrB_Index numEdges = 0;
    GrB_Matrix_nvals(&numEdges, G->A) ;
    numEdges /= 2;
    uint64_t pSwaps = 0; 
    int info = LAGr_SwapEdges(
        G_new, &pSwaps, G, .70, .10, 
        (uint64_t) (numEdges * Q), 891234789234ull, msg) ;
    *pQ = ((double) pSwaps / (double) numEdges);
    return  info;
}
