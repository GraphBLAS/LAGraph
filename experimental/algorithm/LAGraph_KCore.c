//------------------------------------------------------------------------------
// LAGraph_KCore: Single K-core Decomposition Using the GraphBLAS API
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Pranav Konduri, Texas A&M University

//------------------------------------------------------------------------------

// The input is an undirected graph, or a directed graph with a symmetric
// adjacency matrix.  Edge weights are ignored.  On output, decomp(i) = k if
// node i is in the k-core, or empty otherwise.

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&deg) ;               \
    GrB_free (&q) ;                 \
    GrB_free (&delta) ;             \
    GrB_free (&done) ;              \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK                    \
    GrB_free (decomp) ;             \
}

#include "LG_internal.h"

// TODO: revise and add to src
// TODO: need both basic and expert methods; this is mixed
// vanilla OK: no GxB used here

int LAGraph_KCore  // TODO: LAGr_KCore (expert), cache is_symmetric_structure
                   // TODO: cache nself_edges
                   // TODO: cache out degree
(
    // outputs:
    GrB_Vector *decomp,     // kcore decomposition
    // inputs:
    LAGraph_Graph G,        // input graph
    uint64_t k,             // kcore to compute
    char *msg
)
{
    LG_CLEAR_MSG ;

    // declare items
    GrB_Matrix A = NULL;
    GrB_Vector deg = NULL, q = NULL, done = NULL, delta = NULL;

    LG_ASSERT (decomp != NULL, GrB_NULL_POINTER) ;
    (*decomp) = NULL ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->is_symmetric_structure == LAGraph_TRUE))
    {
         // the structure of A is known to be symmetric
        A = G->A ;
    }
    else
    {
        // A is not known to be symmetric
        LG_ASSERT_MSG (false, -1005, "G->A must be symmetric") ;
    }

    // no self edges can be present
    LG_ASSERT_MSG (G->nself_edges == 0, -1004, "G->nself_edges must be zero") ;

    //create work scalars
    GrB_Index n, qnvals, degnvals, maxDeg;
    GRB_TRY (GrB_Matrix_nrows(&n, A)) ;

    //create deg vector using rowdegree property
    LG_TRY (LAGraph_Cached_OutDegree(G, msg)) ;
    GRB_TRY (GrB_Vector_dup(&deg, G->out_degree)) ; //original deg vector is technically 1-core since 0 is omitted
    GRB_TRY (GrB_Vector_nvals(&degnvals, deg)) ;

    //retrieve the max degree level of the graph
    GRB_TRY (GrB_reduce(&maxDeg, GrB_NULL, GrB_MAX_MONOID_INT64, G->out_degree, GrB_NULL)) ;

    //select int type for work vectors and semirings
    GrB_Type int_type  = (maxDeg > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;

    GRB_TRY (GrB_Vector_new(&q, int_type, n));
    GRB_TRY (GrB_Vector_new(&done, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new(&delta, int_type, n)) ;
    //set output to int64
    GRB_TRY (GrB_Vector_new(decomp, int_type, n)) ;

    //change deg vector to int32 if needed
    if (int_type == GrB_INT32)
    {
        // TODO: deg is freed; it should never have been constructed above
        GrB_free (&deg) ;
        GRB_TRY (GrB_Vector_new(&deg, int_type, n)) ;
        GRB_TRY (GrB_assign (deg, G->out_degree, NULL, G->out_degree, GrB_ALL, n, NULL)) ;
    }

    // determine semiring types
    GrB_IndexUnaryOp valueLT = (maxDeg > INT32_MAX) ? GrB_VALUELT_INT64 : GrB_VALUELT_INT32 ;
    GrB_BinaryOp minus_op = (maxDeg > INT32_MAX) ? GrB_MINUS_INT64 : GrB_MINUS_INT32 ;
    GrB_Semiring semiring = (maxDeg > INT32_MAX) ? LAGraph_plus_one_int64 : LAGraph_plus_one_int32 ;

    GRB_TRY (LG_SET_FORMAT_HINT (done, LG_BITMAP + LG_FULL)) ;

    // Creating q
    // get all nodes with degree = level
    GRB_TRY (GrB_select (q, GrB_NULL, GrB_NULL, valueLT, deg, k, GrB_NULL)) ;
    GRB_TRY (GrB_Vector_nvals(&qnvals, q));

    // while q not empty
    int round = 0;
    while (qnvals > 0 && degnvals > 0)
    {
        round++;
        //add anything in q as true into the done list
        //structure to take care of 0-node cases
        GRB_TRY (GrB_assign (done, q, NULL, (bool) true, GrB_ALL, n, GrB_DESC_S)) ;

        // Create delta (the nodes who lost friends, and how many they lost) (push version)
        GRB_TRY (GrB_vxm (delta, GrB_NULL, GrB_NULL, semiring, q, A, GrB_NULL));

        // Create new deg vector (keep anything not in done vector)
        GRB_TRY (GrB_eWiseAdd(deg, done, GrB_NULL, minus_op, deg, delta, GrB_DESC_RSC)) ;

        // Update q, set new nvals
        GRB_TRY (GrB_select (q, GrB_NULL, GrB_NULL, valueLT, deg, k, GrB_NULL)) ;
        GRB_TRY (GrB_Vector_nvals(&qnvals, q)) ;
        GRB_TRY (GrB_Vector_nvals(&degnvals, deg)) ;
    }
    //Assign values of deg to decomp
    GRB_TRY (GrB_assign (*decomp, deg, NULL, k, GrB_ALL, n, GrB_NULL)) ;
    GRB_TRY(GrB_Vector_wait(*decomp, GrB_MATERIALIZE));
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
