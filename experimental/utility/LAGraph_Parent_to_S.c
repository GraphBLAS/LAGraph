//------------------------------------------------------------------------------
// LAGraph_Parent_to_S: Given a dense parent vector for an undirected graph, builds the
// corresponding S matrix needed to coarsen the graph
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2021, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Vidith Madhu, Texas A&M University

//--------

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_Parent_to_S
(
    GrB_Matrix *result,
    GrB_Vector parent,  // dense vector of size n. parent[i] -> representative of node i
    int compress,       // whether to compress the namespace of the coarsened graph
    char *msg
)
{
    GrB_Matrix S ;
    if (compress) {
        
    } else {

    }
    (*result) = S ;
    return (GrB_SUCCESS) ;
}