//------------------------------------------------------------------------------
// LAGraph_MinCut: min cut
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2026 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Darin Peries and Tim Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_MinCut is a GraphBLAS implementation of the push-relabel algorithm
// of Baumstark et al. [1], for computing the minimum cut.
//
// [1] N. Baumstark, G. E. Blelloch, and J. Shun, "Efficient Implementation of
// a Synchronous Parallel Push-Relabel Algorithm." In: Bansal, N., Finocchi, I.
// (eds) Algorithms - ESA 2015. Lecture Notes in Computer Science(), vol 9294.
// Springer, Berlin, Heidelberg.  https://doi.org/10.1007/978-3-662-48350-3 10.

// [2] D. Peries and T. Davis, "A parallel push-relabel maximum flow algorithm
// in LAGraph and GraphBLAS", IEEE HPEC'25, Sept 2025.

// [10] A. V. Goldberg and R. E. Tarjan, “A new approach to the maximum
// flow problem,” in Proc. 18th Annual ACM Symp. Theory of Computing,
// STOC ’86, p. 136–146, ACM, 1986.

// The algorithm takes the residual graph produced by the max flow and runs a
// single breadth-first search from the original source node used in the
// Max Flow. The algorithm returns the set of nodes, S, which are reachable
// from the source, and the set of nodes, S_bar, that are not. Additionally
// The algorithm returns the set of weigthed edges that are included in cut.


// Example: First use the max flow algorithm on a weighted graph with no
// negative edge weights. The optional parameter to generate the R matrix must
// be set. Then take the original graph, the R matrix, and the original
// source node and input them into the parameters of the algorithm.


//------------------------------------------------------------------------------

#include <LAGraph.h>
#include "LG_internal.h"
#include <LAGraph.h>

#undef LG_FREE_ALL
#undef LG_FREE_WORK

#define LG_FREE_WORK				\
{						\
 GrB_free(&S_diag);				\
 GrB_free(&S_bar_diag);				\
 GrB_free(&S_bfs) ;     			\
 G->A = NULL ;  /* do not delete; this is R */  \
 LAGraph_Delete(&G, msg);			\
}

#define LG_FREE_ALL                             \
  { LG_FREE_WORK }

// FIXME: make s input GrB_Scalar (remove t)

int LAGraph_MinCut
(
    // outputs
    GrB_Vector* S,
    GrB_Vector* S_bar,
    GrB_Matrix* cut_set,
    // inputs
    GrB_Matrix R,
    LAGraph_Graph G_origin, // swap G_origin and R
    GrB_Index src,    // src: FIXME
    char *msg
)
{
  //do a bfs from the source to the sink, stop if the frontier is empty 

  LAGraph_Graph G = NULL;
  GrB_Matrix S_diag=NULL, S_bar_diag=NULL ;
  GrB_Vector S_bfs = NULL ;
  GrB_Index n = 0;
  GrB_Matrix_nrows(&n, R);

  LG_TRY(LAGraph_CheckGraph(G_origin, msg));
  LG_ASSERT (S != NULL, GrB_NULL_POINTER) ;
  LG_ASSERT (S_bar != NULL, GrB_NULL_POINTER) ;
  LG_ASSERT (cut_set != NULL, GrB_NULL_POINTER) ;
  LG_ASSERT (src < n, GrB_INVALID_VALUE) ;

  GrB_Matrix A = G_origin->A;
  
  LG_TRY(GrB_Matrix_new(cut_set, GrB_FP64, n, n));
  LG_TRY(GrB_Vector_new(S_bar, GrB_FP64, n));
  LG_TRY(GrB_Vector_new(S, GrB_FP64, n));

  LG_TRY(LAGraph_New(&G, &R, LAGraph_ADJACENCY_DIRECTED, msg));

  //S_bfs is allocated during the bfs
  LG_TRY(LAGr_BreadthFirstSearch(&S_bfs, NULL, G, src, msg));

  LG_TRY(GrB_assign(*S_bar, S_bfs, NULL, 1, GrB_ALL, n, GrB_DESC_SC));
  LG_TRY(GrB_assign(*S, S_bfs, NULL, 1, GrB_ALL, n, GrB_DESC_S));

  GRB_TRY(GrB_Matrix_diag(&S_diag, *S, 0));
  GRB_TRY(GrB_Matrix_diag(&S_bar_diag, *S_bar, 0));

  GRB_TRY(GrB_mxm(*cut_set, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A, S_bar_diag, NULL));
  GRB_TRY(GrB_mxm(*cut_set, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, S_diag, *cut_set, NULL));

  
  LG_FREE_ALL;
  return (GrB_SUCCESS);
}
