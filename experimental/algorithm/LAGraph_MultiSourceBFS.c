//------------------------------------------------------------------------------
// LAGraph_MultiSourceBFS: BFS from several source nodes in parallel
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// TODO: almost ready for src; need a vanilla method

// Takes in a vector of source nodes and finds level and/or parent vectors for each,
// stored together in a matrix

#define LG_FREE_WORK        \
{                           \
    GrB_free (&q) ;         \
}

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&pi) ;        \
    GrB_free (&v) ;         \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_MultiSourceBFS
(
    // outputs:
    GrB_Matrix    *level,
    GrB_Matrix    *parent,
    // inputs:
    const LAGraph_Graph G,
    GrB_Vector      src,
    char          *msg
)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Matrix q = NULL ;           // the current frontier 
    // GrB_Vector w = NULL ;     to compute work remaining, removed since not doing push-pull
    GrB_Matrix pi = NULL ;          // parent matrix
    GrB_Matrix v = NULL ;           // level matrix

    bool compute_level  = (level != NULL) ;
    bool compute_parent = (parent != NULL) ;
    if (compute_level ) (*level ) = NULL ;
    if (compute_parent) (*parent) = NULL ;
    LG_ASSERT_MSG (compute_level || compute_parent, GrB_NULL_POINTER,
        "either level or parent must be non-NULL") ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    //--------------------------------------------------------------------------
    // get the problem size and cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    
    GrB_Index nsrc; // holds the number of source nodes
    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    GRB_TRY (GrB_Vector_size (&nsrc, src)) ;
    for (int64_t s = 0; s < nsrc; s++) 
    {
        GrB_Index currsrc;
        GRB_TRY (GrB_Vector_extractElement (&currsrc, src, s)) ;
        LG_ASSERT_MSG (currsrc < n, GrB_INVALID_INDEX, "invalid source node") ;
    }
    bool very_sparse = (nvals <= n/16) ;

    // determine the semiring type 
    GrB_Type int_type = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    GrB_Semiring semiring ;

    if (compute_parent)
    {
        // use the ANY_SECONDI_INT* semiring: either 32 or 64-bit depending on
        // the # of nodes in the graph.
        semiring = (n > INT32_MAX) ?
            GxB_ANY_SECONDI_INT64 : GxB_ANY_SECONDI_INT32 ;

        // create the parent matrix.  pi(i, j) is the parent id of node j in source i's BFS
        GRB_TRY (GrB_Matrix_new (&pi, int_type, nsrc, n)) ;
        if (!very_sparse)
        {
            GRB_TRY (LG_SET_FORMAT_HINT (pi, LG_BITMAP + LG_FULL)) ;
        }

        // pi (i, src) = src denotes the root of that row's BFS tree
        for (int64_t s = 0; s < nsrc; s++) 
        {
            GrB_Index currsrc;
            GRB_TRY (GrB_Vector_extractElement (&currsrc, src, s)) ;
            // now s contains the row associated with the current source node
            // and currsrc contains the column that should index itself
            GRB_TRY (GrB_Matrix_setElement (pi, currsrc, s, currsrc)) ;
        }

        // create a sparse integer matrix q, and set q(s, src) = src for each row's source
        GRB_TRY (GrB_Matrix_new (&q, int_type, nsrc, n)) ;
        for (int64_t s = 0; s < nsrc; s++) 
        {
            GrB_Index currsrc;
            GRB_TRY (GrB_Vector_extractElement (&currsrc, src, s)) ;
            // now s contains the row associated with the current source node
            // and currsrc contains the column that should index itself
            GRB_TRY (GrB_Matrix_setElement (q, currsrc, s, currsrc)) ;
        }
    }
    else
    {
        // only the level is needed, use the LAGraph_any_one_bool semiring
        semiring = LAGraph_any_one_bool ;

        // create a sparse boolean matrix q, and set q(s, src) = true for the source in each row
        GRB_TRY (GrB_Matrix_new (&q, GrB_BOOL, nsrc, n)) ;
        for (int64_t s = 0; s < nsrc; s++) 
        {
            GrB_Index currsrc;
            GRB_TRY (GrB_Vector_extractElement (&currsrc, src, s)) ;
            // now s contains the row associated with the current source node
            // and currsrc contains the column that should be true
            GRB_TRY (GrB_Matrix_setElement (q, true, s, currsrc)) ;
        }
    }

    if (compute_level)
    {
        // create the level matrix. v(i,j) is the level of node j in source i's BFS
        // v (s, src) = 0 denotes the source node of that row
        GRB_TRY (GrB_Matrix_new (&v, int_type, nsrc, n)) ;
        if (!very_sparse)
        {
            GRB_TRY (LG_SET_FORMAT_HINT (v, LG_BITMAP + LG_FULL)) ;
        }
        for (int64_t s = 0; s < nsrc; s++) 
        {
            GrB_Index currsrc;
            GRB_TRY (GrB_Vector_extractElement (&currsrc, src, s)) ;
            // now s contains the row associated with the current source node
            // and currsrc contains the column that should be 0
            GRB_TRY (GrB_Matrix_setElement (v, 0, s, currsrc)) ;
        }
    }

    GrB_Index nq = nsrc ;      // number of nodes in the current level
    // skipping work remaining computation set-up since we're not doing push-pull

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    // {!mask} is the set of unvisited nodes
    GrB_Matrix mask = (compute_parent) ? pi : v ;

    for (int64_t nvisited = nsrc, k = 1 ; nvisited < n*nsrc ; nvisited += nq, k++)
    {

        //----------------------------------------------------------------------
        // q = frontier at the kth level of the BFS
        //----------------------------------------------------------------------

        // mask is pi if computing parent, v if computing just level
        
        // push (saxpy-based mxm):  q'{!mask} = q'*A 
        GRB_TRY (GrB_mxm (q, mask, NULL, semiring, q, A, GrB_DESC_RSC)) ;
        
        //----------------------------------------------------------------------
        // done if q is empty
        //----------------------------------------------------------------------

        GRB_TRY (GrB_Matrix_nvals (&nq, q)) ;
        if (nq == 0)
        {
            break ;
        }

        //----------------------------------------------------------------------
        // assign parents/levels
        //----------------------------------------------------------------------

        if (compute_parent)
        {
            // q(s, i) currently contains the parent id of node i in tree s.
            // pi{q} = q
            GRB_TRY (GrB_assign (pi, q, NULL, q, GrB_ALL, nsrc, GrB_ALL, n, GrB_DESC_S)) ;
        }
        if (compute_level)
        {
            // v{q} = k, the kth level of the BFS
            GRB_TRY (GrB_assign (v, q, NULL, k, GrB_ALL, nsrc, GrB_ALL, n, GrB_DESC_S)) ;
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    if (compute_parent) (*parent) = pi ;
    if (compute_level ) (*level ) = v ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
