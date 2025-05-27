//------------------------------------------------------------------------------
// LAGr_PageRankGX: pagerank for the LDBC Graphalytics (GX) benchmark
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis and Mohsen Aznaveh, Texas A&M University;
// Originally contributed by Gabor Szarnyas and Balint Hegyi, Budapest
// University of Technology and Economics (with accented characters: G\'{a}bor
// Sz\'{a}rnyas and B\'{a}lint Hegyi, using LaTeX syntax).

//------------------------------------------------------------------------------

// This is an Advanced algorithm (G->AT and G->out_degree are required).

// PageRank for the LDBC Graphalytics (GX) benchmark.
// Do not use in production.

// This algorithm follows the specification given in the LDBC Graphalytics
// benchmark, see https://arxiv.org/pdf/2011.15028.pdf

// The G->AT and G->out_degree cached properties must be defined for this
// method.  If G is undirected or G->A is known to have a symmetric structure,
// then G->A is used instead of G->AT, however.

#define LG_FREE_WORK                \
{                                   \
    GrB_free (&d1) ;                \
    GrB_free (&d) ;                 \
    GrB_free (&t) ;                 \
    GrB_free (&w) ;                 \
    GrB_free (&non_sink_mask) ;     \
    GrB_free (&sink_vec) ;          \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (&r) ;                 \
}

#include "LG_internal.h"

int LAGr_PageRankGX
(
    // output:
    GrB_Vector *centrality, // centrality(i): GX-style pagerank of node i
    int *iters,             // number of iterations taken
    // input:
    const LAGraph_Graph G,  // input graph
    float damping,          // damping factor (typically 0.85)
    int itermax,            // maximum number of iterations (typically 100)
    char *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector r = NULL, d = NULL, t = NULL, w = NULL, d1 = NULL ;
    GrB_Vector non_sink_mask = NULL, sink_vec = NULL;

    LG_ASSERT (centrality != NULL, GrB_NULL_POINTER) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GrB_Matrix AT ;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        AT = G->A ;
    }
    else
    {
        // A and A' differ
        AT = G->AT ;
        LG_ASSERT_MSG (AT != NULL,
            LAGRAPH_NOT_CACHED, "G->AT is required") ;
    }
    GrB_Vector d_out = G->out_degree ;
    LG_ASSERT_MSG (d_out != NULL,
        LAGRAPH_NOT_CACHED, "G->out_degree is required") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n ;
    (*centrality) = NULL ;
    GRB_TRY (GrB_Matrix_nrows (&n, AT)) ;

    const float scaled_damping = (1 - damping) / n ;
    const float teleport = scaled_damping ; // teleport = (1 - damping) / n

    // Determine vector of non-sink vertices:
    // These vertices are the ones which have outgoing edges. In subsequent
    // operations, this mask can be negated to select sink vertices.
    GRB_TRY (GrB_Vector_new (&non_sink_mask, GrB_BOOL, n)) ;
    GRB_TRY (GrB_reduce (non_sink_mask, NULL, NULL, GrB_LOR_MONOID_BOOL,
        G->A, NULL)) ;

    // Initialize vector for collecting the sink values
    GRB_TRY (GrB_Vector_new (&sink_vec, GrB_FP64, n)) ;

    // r = 1 / n
    GRB_TRY (GrB_Vector_new (&t, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&r, GrB_FP64, n)) ;
    GRB_TRY (GrB_Vector_new (&w, GrB_FP64, n)) ;
    GRB_TRY (GrB_assign (r, NULL, NULL, (float) (1.0 / n), GrB_ALL, n, NULL)) ;

    // prescale with damping factor, so it isn't done each iteration
    // d = d_out / damping ;
    GRB_TRY (GrB_Vector_new (&d, GrB_FP64, n)) ;
    GRB_TRY (GrB_apply (d, NULL, NULL, GrB_DIV_FP64, d_out, damping, NULL)) ;

    // d1 = 1 / damping
    float dmin = 1.0 / damping ;
    GRB_TRY (GrB_Vector_new (&d1, GrB_FP64, n)) ;
    GRB_TRY (GrB_assign (d1, NULL, NULL, dmin, GrB_ALL, n, NULL)) ;
    // d = max (d1, d)
    GRB_TRY (GrB_eWiseAdd (d, NULL, NULL, GrB_MAX_FP64, d1, d, NULL)) ;
    GrB_free (&d1) ;

    //--------------------------------------------------------------------------
    // pagerank iterations
    //--------------------------------------------------------------------------

    for ((*iters) = 0 ; (*iters) < itermax ; (*iters)++)
    {
        // swap t and r ; now t is the old score
        GrB_Vector temp = t ; t = r ; r = temp ;

        // sink value calculation
        // Extract all the sink PR entries from the previous result
        GRB_TRY (GrB_extract (sink_vec, non_sink_mask, NULL, t, GrB_ALL,
            n, GrB_DESC_SC)) ;

        // Sum the previous PR values of sink vertices together
        double sink_value;
        GRB_TRY (GrB_reduce (&sink_value, NULL, GrB_PLUS_MONOID_FP64,
            sink_vec, NULL )) ;

        // Multiply by damping factor and 1 / |V|
        sink_value *= (damping / n);

        // w = t ./ d
        GRB_TRY (GrB_eWiseMult (w, NULL, NULL, GrB_DIV_FP64, t, d, NULL)) ;
        // r = teleport + redistributed from sinks
        GRB_TRY (GrB_assign (r, NULL, NULL, teleport + sink_value, GrB_ALL,
            n, NULL)) ;
        // r += A'*w
        GRB_TRY (GrB_mxv (r, NULL, GrB_PLUS_FP64, LAGraph_plus_second_fp64,
            AT, w, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*centrality) = r ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
