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
    LAGraph_Graph *G_new, //The adjacency matrix of G with edges randomly swapped
    // input: not modified
    LAGraph_Graph G,
    GrB_Index Q, // Swaps per edge
    char *msg
)
{
    GrB_Index numSwaps = 0;
    GrB_Matrix_nvals(&numSwaps, G->A) ;
    numSwaps /= 2;
    numSwaps *= Q;
    return LAGr_SwapEdges(G_new, G, .70, .10, numSwaps, 891234789234ull, msg) ;
}
