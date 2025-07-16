//------------------------------------------------------------------------------
// LAGraph_EstimateDiameter: Graph Diameter Estimation
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Alexandra Goff

//------------------------------------------------------------------------------

// TODO: almost ready for src; need to handle GxB

// Takes in a graph and estimates the diameter 
// and optionally also finds pseudo-peripheral nodes of the graph

// Outputs: 
// Diameter returns the estimated diameter of the graph
// If not set to NULL, peripheral will be a vector with n elements
// index i of peripheral is the estimated diameter if it's a pseudo-peripheral node or nothing if not

// Inputs:
// G is the graph to be analyzed
// maxSrcs limits the number of sources used each cycle
// maxLoops limits the number of times the core loop will run if a stable diameter isn't found
// msg is a buffer for error messages

#define LG_FREE_WORK        \
{                           \
    GrB_free (&ecc) ;       \
    GrB_free (&candidateSrcs) ;       \
    GrB_free (&srcs) ;      \
    GrB_free (&level) ;     \
    GrB_free (&Mod) ;       \
    LAGraph_Free((void**) &sourceIndicies, NULL);   \
    LAGraph_Free((void**) &sourceValues, NULL);     \
}

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&peri) ;      \
}

#include "LG_internal.h"
#include "LAGraphX.h"

void LG_ED_mod32 (int32_t *z, const int32_t *x, const int32_t *y) ;
void LG_ED_mod32 (int32_t *z, const int32_t *x, const int32_t *y)
{
    /* make sure x is positive */
    int32_t t = ((*x) > 0) ? (*x) : -(*x) ;
    (*z) = t % (*y) ;
}

#define MOD32_DEFN                                                  \
"void LG_ED_mod32 (int32_t *z, const int32_t *x, const int32_t *y)    \n" \
"{                                                              \n" \
"    /* make sure x is positive */                              \n" \
"    int32_t t = ((*x) > 0) ? (*x) : -(*x) ;                    \n" \
"    (*z) = t % (*y) ;                                          \n" \
"}"

void LG_ED_mod64 (int64_t *z, const int64_t *x, const int64_t *y) ;
void LG_ED_mod64 (int64_t *z, const int64_t *x, const int64_t *y)
{
    /* make sure x is positive */
    int64_t t = ((*x) > 0) ? (*x) : -(*x) ;
    (*z) = t % (*y) ;
}

#define MOD64_DEFN                                                  \
"void LG_ED_mod64 (int64_t *z, const int64_t *x, const int64_t *y)    \n" \
"{                                                              \n" \
"    /* make sure x is positive */                              \n" \
"    int64_t t = ((*x) > 0) ? (*x) : -(*x) ;                    \n" \
"    (*z) = t % (*y) ;                                          \n" \
"}"

int LAGraph_EstimateDiameter
(
    // outputs:
    GrB_Index    *diameter,
    GrB_Vector    *peripheral,
    // inputs:
    const LAGraph_Graph G,
    GrB_Index    maxSrcs,
    GrB_Index    maxLoops,
    uint64_t     seed,          // seed for randomization
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
    GrB_Index d = 0 ;              // current diameter
    GrB_Index lastd = 0 ;          // previous diameter
    GrB_Vector srcs = NULL ;        // list of current sources
    GrB_Index nsrcs ;               // number of current sources
    GrB_Matrix level = NULL ;       // matrix for msbfs to put level info in
    GrB_Vector candidateSrcs = NULL ; // work vector for getting sources for the next iteration of the loop
    GrB_BinaryOp Mod = NULL ;

    GrB_Index *sourceIndicies = NULL ;
    int64_t *sourceValues = NULL ; // just need this so extractTuples will run

    bool compute_periphery  = (peripheral != NULL) ;
    if (compute_periphery ) (*peripheral) = NULL ;
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

    //--------------------------------------------------------------------------
    // set up the first maxSrcs random nodes
    //--------------------------------------------------------------------------

    // currently just doing the first maxSrcs, consider different randomization
    // check maxSrcs < n
    if (maxSrcs > n)
    {
        nsrcs = n;
    }
    else
    {
        nsrcs = maxSrcs;
    }
    GRB_TRY (GrB_Vector_new (&srcs, int_type, nsrcs)) ;
    // srcs (0:nsrcs-1) = 0
    GRB_TRY (GrB_assign (srcs, NULL, NULL, 0, GrB_ALL, nsrcs, NULL)) ;
    if (nsrcs == n)
    {
        // srcs = 0:n-1
        GrB_IndexUnaryOp op =
            (n > INT32_MAX) ?  GrB_ROWINDEX_INT64 : GrB_ROWINDEX_INT32 ;
        GRB_TRY (GrB_apply (srcs, NULL, NULL, op, srcs, 0, NULL)) ;
    }
    else
    {
        // srcs = randomized, still of size nsrcs-1, values in range 0 to UINT64_MAX
        // TODO: if the graph is very sparse, select nodes at random with at
        // least one out-going edge.
        LAGRAPH_TRY (LAGraph_Random_Seed (srcs, seed, msg)) ;
        GRB_TRY (GxB_BinaryOp_new (&Mod,
            (n > INT32_MAX) ? ((GxB_binary_function)LG_ED_mod64) :
                              ((GxB_binary_function)LG_ED_mod32),
            int_type, int_type, int_type,
            (n > INT32_MAX) ? "LG_ED_mod64" : "LG_ED_mod32",
            (n > INT32_MAX) ? MOD64_DEFN : MOD32_DEFN)) ;
        GRB_TRY (GrB_apply (srcs, NULL, NULL, Mod, srcs, n, NULL)) ;
        GrB_free (&Mod) ;
    }

    //--------------------------------------------------------------------------
    // core loop, run until current and previous diameters match or reach given limit
    //--------------------------------------------------------------------------

    GrB_Monoid max =
        (n > INT32_MAX) ?  GrB_MAX_MONOID_INT64 : GrB_MAX_MONOID_INT32 ;
    GrB_IndexUnaryOp eqOp =
        (n > INT32_MAX) ?  GrB_VALUEEQ_INT64 : GrB_VALUEEQ_INT32 ;
    bool incSrcs = false;
    for (int64_t i = 0; i < maxLoops; i++)
    {
        // save previous diameter
        lastd = d;

        // get new diameter 
        GrB_free (&level) ;
        LG_TRY (LAGraph_MultiSourceBFS(&level, NULL, G, srcs, msg)) ;
        GrB_free (&ecc) ;
        GRB_TRY (GrB_Vector_new (&ecc, int_type, n)) ;
        GRB_TRY (GrB_reduce(ecc, NULL, NULL, max, level, GrB_DESC_T0)) ;
        GRB_TRY (GrB_reduce(&d, NULL, max, ecc, GrB_NULL)) ;
        GrB_free (&level) ;

        // check if done
        if (d == lastd){
            incSrcs = true;
            break;
        }

        // now with fewer for loops
        // said in last discussion: remaining for loop fine because of batch processing?

        // set up source list for next round
        // get the number of peripheral nodes 
        GrB_Index nperi = 0;
        GRB_TRY (GrB_Vector_new (&candidateSrcs, int_type, n)) ;
        GRB_TRY (GrB_select(candidateSrcs, NULL, NULL, eqOp, ecc, d, NULL)) ;
        GRB_TRY (GrB_Vector_nvals(&nperi, candidateSrcs)) ;
        

        // select the number of sources
        if (nperi > maxSrcs) {
            nsrcs = maxSrcs;
        } else {
            nsrcs = nperi;
        }

        // choose sources
        GrB_free (&srcs) ;
        GRB_TRY (GrB_Vector_new (&srcs, int_type, nsrcs)) ;

        GRB_TRY (LAGraph_Calloc((void **) &sourceIndicies, nperi, sizeof(GrB_Index), msg)) ;
        GRB_TRY (LAGraph_Calloc((void **) &sourceValues, nperi, sizeof(int64_t), msg)) ;
        GRB_TRY (GrB_Vector_extractTuples(sourceIndicies, sourceValues, &nperi, candidateSrcs)) ;
        for (int64_t j = 0; j < nsrcs; j++) {
            GRB_TRY (GrB_Vector_setElement (srcs, sourceIndicies[j], j)) ;
        }
        LAGraph_Free((void**) &sourceIndicies, NULL);
        LAGraph_Free((void**) &sourceValues, NULL);
        GrB_free(&candidateSrcs) ;
    }

    //--------------------------------------------------------------------------
    // after loop, set up peripheral nodes if needed
    //--------------------------------------------------------------------------
    if (compute_periphery) {
        GRB_TRY (GrB_Vector_new (&peri, int_type, n)) ;

        GRB_TRY (GrB_select(peri, NULL, NULL, eqOp, ecc, d, NULL)) ;

        if (incSrcs) {
            for (int64_t i = 0; i < nsrcs; i++) {
                GrB_Index currsrc;
                GRB_TRY(GrB_Vector_extractElement(&currsrc, srcs, i));
                GRB_TRY (GrB_Vector_setElement (peri, d, currsrc)) ;  
            }
        }
       

    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    if (compute_periphery) (*peripheral) = peri ;
    (*diameter ) = d ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
