//------------------------------------------------------------------------------
// LAGraph_MaximumMatching: maximum matching between nodes of disjoint sets
//                          in bipartite graphs
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Christina Koutsou, Aristotle University of Thessaloniki

// add paper

//------------------------------------------------------------------------------

// add explanation of paper

// #define LG_FREE_WORK \
//     {
// }

// #define LG_FREE_ALL \
//     {               \
//         // LG_FREE_WORK ;                  \
//     // GrB_free (centrality) ;
// }

#include "LG_internal.h"

//------------------------------------------------------------------------------
// the Vertex tuple: (parentC, rootC)
//------------------------------------------------------------------------------

typedef struct
{
    uint64_t parentC;
    uint64_t rootC;
} vertex;

// repeat the typedef as a string, to give to GraphBLAS
#define VERTEX_DEFN      \
    "typedef struct "    \
    "{ "                 \
    "uint64_t parentC; " \
    "uint64_t rootC; "   \
    "} "                 \
    "vertex; "

void *initFrontier(vertex *z, vertex *x, GrB_Index i, GrB_Index j, const void *y)
{

    x->parentC = i;
    x->rootC = i;
    z = x;
}

#define INIT_FRONTIER_DEFN                                                      \
    "vertex *initFrontier(vertex *x, GrB_Index i, GrB_Index j, const void *y) " \
    "{ "                                                                        \
    "x->parentC = i; "                                                          \
    "x->rootC = i; "                                                            \
    "return x; "                                                                \
    "} "                                                                        \
    "vertex;"

int LAGraph_MaximumMatching(
    // output/input:
    GrB_Vector *mateC, // mateC(j) = i : Column j of C subset is matched to row i of R subset
    // input:
    LAGraph_Graph G, // input graph
    char *msg)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix A = G->A;
    uint64_t ncols = 0;
    GRB_TRY(GrB_Matrix_ncols(&ncols, A));

    uint64_t nrows = 0;
    GRB_TRY(GrB_Matrix_nrows(&nrows, A));

    GrB_Vector pathC = NULL; // make this bitmap, if dense I would have to give all the entries and make the matrix 1-based
    GRB_TRY(GrB_Vector_new(&pathC, GrB_UINT64, ncols));

    GrB_Vector parentsR = NULL; // parents of row nodes that are reachable from paths of the initial column frontier
    GRB_TRY(GrB_Vector_new(&parentsR, GrB_UINT64, nrows));

    GrB_Type Vertex = NULL;
    GRB_TRY(GxB_Type_new(&Vertex, sizeof(vertex), "vertex", VERTEX_DEFN));

    GrB_Vector frontierC = NULL;
    GRB_TRY(GrB_Vector_new(&frontierC, Vertex, ncols));
    GrB_Vector frontierR = NULL;
    GRB_TRY(GrB_Vector_new(&frontierR, Vertex, nrows));

    GrB_IndexUnaryOp initFrontierOp = NULL;
    GRB_TRY(GxB_IndexUnaryOp_new(&initFrontierOp, (GxB_index_unary_function)initFrontier, Vertex, Vertex, NULL, "initFrontier", INIT_FRONTIER_DEFN));

    do
    {
        GRB_TRY(GrB_Vector_clear(pathC));
        GRB_TRY(GrB_Vector_clear(parentsR));

        // for every col j not matched, assign f(j) = VERTEX(j,j)
        GRB_TRY(GrB_Vector_apply_IndexOp_UDT(frontierC, mateC, NULL, initFrontierOp, frontierC, NULL, GrB_DESC_RSC));

    } while (pathC != NULL); // only in the first and last iteration should this condition be false
}