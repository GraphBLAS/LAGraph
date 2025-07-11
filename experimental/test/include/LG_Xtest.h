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

int LG_check_kcore
(
    // outputs:
    GrB_Vector *decomp,     // kcore decomposition
    uint64_t *kmax,         // max kcore- if kfinal == -1, kmax = -1
    // inputs
    LAGraph_Graph G,        // input graph
    int64_t kfinal,         // max k to check for graph.
    char *msg
) ;

int LG_check_kcore_decompose
(
    // outputs:
    GrB_Matrix *D,              // kcore decomposition
    // inputs:
    LAGraph_Graph G,            // input graph
    GrB_Vector decomp,
    uint64_t k,
    char *msg
) ;

int LG_check_lcc
(
     // outputs:
     GrB_Vector *coefficients,     // the local clustering coefficients
     // inputs
     LAGraph_Graph G,        // input graph
     char *msg
) ;

int LG_check_coarsen
(
    // outputs:
    GrB_Matrix *coarsened,    // coarsened adjacency
    // inputs:
    GrB_Matrix A,               // input adjacency (for the purposes of testing, is FP64)
    GrB_Vector parent,          // parent mapping. Must not be NULL.
    GrB_Vector newlabel,       // new labels of nodes, used to populate resulting adjacency matrix, can be NULL if preserve_mapping = 1, else must be a valid result
    GrB_Vector inv_newlabel,   // inverse of newlabel, can be NULL if preserve_mapping = 1, else must be a valid result
    int preserve_mapping,       // whether to preserve the original namespace of nodes
    int combine_weights,        // whether to combine the weights of edges that collapse together
    char *msg
) ;

int LG_check_coloring
(
    // inputs
    LAGraph_Graph G,
    GrB_Vector C,
    char *msg
) ;

int LG_check_edgeBetweennessCentrality
(
    // output
    GrB_Matrix *C,       // centrality matrix
    // input
    LAGraph_Graph G,
    GrB_Vector sources,  // source vertices to compute shortest paths (if NULL or empty, use all vertices)
    char *msg
) ;

int LG_check_rcc
(
    GrB_Vector *rich_club_coefficents, //output
    LAGraph_Graph G, //input graph
    char *msg
) ;

int LG_check_argminmax
(
    // output
    GrB_Vector *x_result,       // min/max value in each row/col of A
    GrB_Vector *p_result,       // index of min/max value in each row/col of A
    // input
    GrB_Matrix A,
    int dim,                    // dim=1: cols of A, dim=2: rows of A
    bool is_min,
    char *msg
) ;

int LG_check_flow
(
   const GrB_Matrix flow_mtx,
   char* msg
) ;

#endif
