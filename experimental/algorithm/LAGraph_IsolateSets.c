//------------------------------------------------------------------------------
// LAGraph_IsolateSets.c: Returns an Isolate Set
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2024 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Olumayowa Olowomeye, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#undef LG_FREE_ALL
#define LG_FREE_ALL               \
    {                             \
        GrB_free(&score);         \
        GrB_free(&scoreA);        \
        GrB_free(&neighbor_max);  \
        GrB_free(&new_members);   \
        GrB_free(&new_membersA);  \
        GrB_free(&new_neighbors); \
        GrB_free(&candidates);    \
        GrB_free(&empty);         \
        GrB_free(&Seed);          \
        GrB_free(&degree);        \
    }
#define DEBUG 0
#if DEBUG
#define check() printf("here")
#define dbg(x) \
    if (DEBUG) \
    GxB_print(x, 5)
#define err(x, info)                                    \
    if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) \
    {                                                   \
        char **err;                                     \
        GrB_error(err, x);                              \
        printf("\ninfo: %d error: %s\n", info, err);    \
    }
#else
#define check()
#define dbg(x)
#define err(x, info)
#endif

int LAGraph_IsolateSet(
    // output
    GrB_Vector *isolate_set,
    // input
    GrB_Matrix A,
    GrB_Vector ignore_node,
    uint64_t seed,
    char *msg)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10_2
    LG_CLEAR_MSG;

    char MATRIX_TYPE[LAGRAPH_MSG_LEN];
    // GrB_set(GrB_GLOBAL, DEBUG, GxB_BURBLE);
    GrB_Vector iset = NULL; // independent set (output vector)
    GrB_Vector score = NULL;
    GrB_Vector scoreA = NULL; // random score for each node
    // GrB_Vector scoreAA = NULL;
    GrB_Vector neighbor_max = NULL; // value of max neighbor score
    GrB_Vector new_members = NULL;  // set of new members to add to iset
    GrB_Vector new_membersA = NULL;
    GrB_Vector new_neighbors = NULL; // new neighbors to new iset members
    GrB_Vector candidates = NULL;    // candidate nodes
    GrB_Vector empty = NULL;         // an empty vector
    GrB_Vector Seed = NULL;          // random number seed vector
    GrB_Vector degree = NULL;        // (float) G->out_degree
    // GrB_Matrix A ;                      // G->A, the adjacency matrix
    GrB_Index n; // # of nodes
    // printf("in Isolate set algorithm");
    // LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT(isolate_set != NULL, GrB_NULL_POINTER);
    // A = G->A;
    dbg(A);

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Vector_new(&iset, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&neighbor_max, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&degree, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&new_members, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&new_neighbors, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&new_membersA, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&candidates, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&empty, GrB_BOOL, n));
    GRB_TRY(GrB_Vector_new(&Seed, GrB_UINT64, n));
    GRB_TRY(GrB_Vector_new(&score, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&scoreA, GrB_FP64, n));

    // rand
    //  seed = 6247;
    //  printf("%ld",seed);
    GRB_TRY(GrB_Matrix_reduce_Monoid(degree, NULL, NULL, GrB_PLUS_MONOID_FP64, A, NULL));
    dbg(degree);

    GrB_Index ncandidates;
    if (ignore_node == NULL)
    {
        GRB_TRY(GrB_assign(candidates, NULL, NULL, (bool)true, GrB_ALL, n, NULL));
        GRB_TRY(GrB_assign(Seed, NULL, NULL, 1, GrB_ALL, n, NULL));
    }
    else
    {
        GRB_TRY(GrB_assign(candidates, ignore_node, NULL, (bool)true, GrB_ALL, n, GrB_DESC_C));
        GRB_TRY(GrB_assign(Seed, candidates, NULL, 1, GrB_ALL, n, GrB_DESC_S));
    }
    dbg(candidates);
    GRB_TRY(LAGraph_Random_Seed(Seed, seed, msg));
    dbg(Seed);

    GRB_TRY(GrB_Vector_nvals(&ncandidates, candidates));

    GRB_TRY(GrB_assign(score, NULL, NULL, Seed, GrB_ALL, n, NULL));
    GRB_TRY(GrB_eWiseMult(score, NULL, NULL, GrB_DIV_FP64, score, degree, NULL));
    dbg(score);

    dbg(candidates);
    GRB_TRY(GrB_vxm(scoreA, candidates, NULL, GrB_MAX_FIRST_SEMIRING_FP64, score, A, GrB_DESC_RS));
    dbg(scoreA);
    GRB_TRY(GrB_vxm(neighbor_max, candidates, NULL, GrB_MAX_FIRST_SEMIRING_FP64, scoreA, A, GrB_DESC_RS));
    GRB_TRY(GrB_vxm(neighbor_max, candidates, NULL, GrB_MAX_FIRST_SEMIRING_FP64, neighbor_max, A, GrB_DESC_RS));
    dbg(neighbor_max);
    dbg(score);

    GRB_TRY(GrB_eWiseAdd(new_members, NULL, NULL, GrB_GE_FP64, score, neighbor_max, NULL));
    dbg(new_members);
    GRB_TRY(GrB_select(new_members, NULL, NULL, GrB_VALUEEQ_BOOL,
                       new_members, (bool)true, NULL));
    dbg(new_members);
    GRB_TRY(GrB_assign(iset, new_members, NULL, true, GrB_ALL, n, NULL));
    (*isolate_set) = iset;

    // printf("done iset");
    iset = NULL;
    dbg(*isolate_set);
    LG_FREE_ALL;
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
    return (GrB_SUCCESS) ;
}

// #undef LG_FREE_WORK
#undef LG_FREE_ALL
#define LG_FREE_ALL

int LAGraph_IsolateSets(
    GrB_Matrix *IsolateSets, // Output: k x n Boolean matrix
    // LAGraph_Graph G,         // Input: graph
    GrB_Matrix A,
    // GrB_Vector ignore_nodes,
    uint64_t seed, // Input: RNG seed
    char *msg      // Error message buffer
)
{
#if LG_SUITESPARSE_GRAPHBLAS_V10_2
    LG_CLEAR_MSG;
    // LG_TRY(LAGraph_CheckGraph(G, msg));
    LG_ASSERT(IsolateSets != NULL, GrB_NULL_POINTER);

    // GrB_Matrix A = G->A;
    GrB_Index n;
    GrB_Vector ignore_nodes = NULL;

    GRB_TRY(GrB_Matrix_nrows(&n, A));
    GRB_TRY(GrB_Vector_new(&ignore_nodes, GrB_BOOL, n));
    GRB_TRY(GrB_assign(ignore_nodes, NULL, NULL, (bool)false, GrB_ALL, n, NULL));

    GrB_Index max_k = n; // Max possible number of isolate sets is <= n
    GrB_Matrix result = NULL;
    GRB_TRY(GrB_Matrix_new(&result, GrB_BOOL, max_k, n));

    GrB_Vector iset = NULL; // start NULL -- LAGraph_IsolateSet will allocate
    GrB_Index k = 0;
    GrB_Index vals_res = 0;

    while (true)
    {
        // LAGraph_IsolateSet allocates a new vector and stores it into 'iset'
        GRB_TRY(LAGraph_IsolateSet(&iset, A, ignore_nodes, seed, msg));

        // if 'iset' is NULL or empty, break
        GRB_TRY(GrB_Vector_nvals(&vals_res, iset));
        if (vals_res == 0)
        {
            GrB_free(&iset); // free the empty iset returned
            break;
        }

        // mark ignored nodes (update ignore_nodes) BEFORE copying into matrix if desired
        GRB_TRY(GrB_Vector_eWiseAdd_BinaryOp(ignore_nodes, NULL, NULL, GrB_LOR, ignore_nodes, iset, NULL));

        // copy the iset vector into the result matrix row k
        GRB_TRY(GxB_Row_assign_Vector(result, NULL, NULL, iset, k, NULL, NULL));

        // free the iset after copying into the matrix to avoid leaking
        GrB_free(&iset);
        iset = NULL;

        k++;
    }
    GRB_TRY(GrB_Matrix_resize(result, k, n));
    *IsolateSets = result;
    GrB_free(&ignore_nodes);
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED);
#endif
}
