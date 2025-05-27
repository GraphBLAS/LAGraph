//------------------------------------------------------------------------------
// LAGraph_ExactDiameter: Graph Diameter Computation
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// TODO: almost ready for src; says it's not vanilla but should be OK

// Takes in a graph and finds the diameter 
// and optionally also the peripheral nodes of the graph

// Outputs: 
// Diameter returns the diameter of the graph
// If not set to NULL, peripheral will be a vector with n elements
// index i of peripheral is the diameter if it's a peripheral node or nothing if not
// If not set to NULL, eccentricity will be a vector with the eccentricity of each 
// node in the graph

// Inputs:
// G is the graph to be analyzed
// k is the number of nodes in each batch of the BFS, 
// higher k value allows for more parallelization at the cost of more space used
// msg is a buffer for error messages

#define LG_FREE_WORK        \
{                           \
    GrB_free (&srcEcc) ;    \
    GrB_free (&ecc) ;       \
    GrB_free (&level) ;     \
    GrB_free (&srcs) ;      \
    LAGraph_Free((void**) &sources, NULL);  \
}

#define LG_FREE_ALL             \
{                               \
    LG_FREE_WORK ;              \
    GrB_free (eccentricity) ;   \
    GrB_free (&peri) ;          \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_ExactDiameter
(
    // outputs:
    GrB_Index    *diameter,
    GrB_Vector    *peripheral,
    GrB_Vector    *eccentricity,
    // inputs:
    const LAGraph_Graph G,
    GrB_Index      k,
    char          *msg
)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector ecc = NULL ;           // the eccentricity of the nodes
    GrB_Vector peri = NULL ;          // vector to store peripheral node status in
    GrB_Vector srcEcc = NULL ;        // work vector to store each iteration's eccentricity in temporarily
    GrB_Matrix level = NULL;          // work matrix for storing msbfs level info
    GrB_Index d ;                     // diameter
    GrB_Vector srcs = NULL ;
    GrB_Index *sources = NULL ;

    bool compute_periphery  = (peripheral != NULL) ;
    if (compute_periphery ) (*peripheral) = NULL ;
    bool compute_eccentricity  = (eccentricity != NULL) ;
    if (compute_eccentricity ) (*eccentricity) = NULL ;
    bool compute_diameter  = (diameter != NULL) ;
    LG_ASSERT_MSG (compute_diameter, GrB_NULL_POINTER,
        "Diameter destination must be non-NULL") ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    //--------------------------------------------------------------------------
    // get the problem size and cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    
    GrB_Index n;        // number of nodes in the graph
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;

    GrB_Type int_type = (n > INT32_MAX) ? GrB_INT64 : GrB_INT32 ;
    GRB_TRY (GrB_Vector_new (&ecc, int_type, n)) ;

    //--------------------------------------------------------------------------
    // get eccentricity, k nodes at a time
    //--------------------------------------------------------------------------

    int64_t setStart = 0;
    GrB_Monoid max = (n > INT32_MAX) ?
            GrB_MAX_MONOID_INT64 : GrB_MAX_MONOID_INT32 ;
    while (setStart < n){
        // set up the sources for this iteration
        int64_t nsrcs;
        if ((setStart + k) <= n){
            nsrcs = k;
        } else {
            nsrcs = n - setStart;
        }
        GRB_TRY (GrB_Vector_new (&srcs, int_type, nsrcs)) ;
        GRB_TRY (LAGraph_Calloc((void **) &sources, nsrcs, sizeof(GrB_Index), msg)) ;
        for (int64_t i = 0; i < nsrcs; i++){
            GRB_TRY (GrB_Vector_setElement (srcs, setStart+i, i)) ;
            sources[i] = setStart+i;
        }

        // run bfs to get level matrix for the sources
        GrB_free (&level) ;
        LAGRAPH_TRY (LAGraph_MultiSourceBFS (&level, NULL, G, srcs, msg)) ;
        GrB_free (&srcs) ;

        // populate setStart to setStart+nsrcs of ecc with the max level for each src
        GRB_TRY (GrB_Vector_new (&srcEcc, int_type, nsrcs)) ;
        GRB_TRY (GrB_reduce(srcEcc, NULL, NULL, max, level, GrB_NULL)) ;
        LAGRAPH_TRY (GrB_assign(ecc, NULL, NULL, srcEcc, sources,  nsrcs, GrB_NULL)) ;
        GrB_free (&level) ;
        GrB_free (&srcEcc) ;

        LAGraph_Free((void**) &sources, NULL);

        // adjust setStart for next iteration
        setStart = setStart + nsrcs;
    }


    //--------------------------------------------------------------------------
    // determine diameter from eccentricity list
    //--------------------------------------------------------------------------
    
    GRB_TRY (GrB_reduce(&d, NULL, max, ecc, GrB_NULL)) ;

    //--------------------------------------------------------------------------
    // get peripheral nodes, if applicable
    //--------------------------------------------------------------------------

    if (compute_periphery){
        GrB_IndexUnaryOp eqOp = (n > INT32_MAX) ?
            GrB_VALUEEQ_INT64 : GrB_VALUEEQ_INT32 ;
        GRB_TRY (GrB_Vector_new (&peri, int_type, n)) ;
        GRB_TRY (GrB_select(peri, NULL, NULL, eqOp, ecc, d, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    if (compute_periphery)
    {
        (*peripheral) = peri ;
        peri = NULL ;
    }
    if (compute_eccentricity)
    {
        (*eccentricity) = ecc ;
        ecc = NULL; // makes sure eccentricity vector doesn't get freed if the user wants it
    } 
    (*diameter) = d;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
