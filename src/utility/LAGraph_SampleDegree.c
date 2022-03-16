//------------------------------------------------------------------------------
// LAGraph_SampleDegree: sample the degree median and mean
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

// LAGraph_SampleDegree computes estimates of the mean and median of the
// row or column degree of a graph.

#define LG_FREE_ALL LAGraph_Free ((void **) &samples, NULL) ;

#include "LG_internal.h"

int LAGraph_SampleDegree
(
    // output:
    double *sample_mean,    // sampled mean degree
    double *sample_median,  // sampled median degree
    // input:
    const LAGraph_Graph G,  // graph of n nodes
    bool byrow,             // if true, sample G->rowdegree, else G->coldegree
    int64_t nsamples,       // number of samples
    uint64_t seed,          // random number seed
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    int64_t *samples = NULL ;
    LG_ASSERT (sample_mean != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (sample_median != NULL, GrB_NULL_POINTER) ;
    nsamples = LAGRAPH_MAX (nsamples, 1) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    GrB_Vector Degree ;

    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGraph_ADJACENCY_DIRECTED &&
        G->structure_is_symmetric == LAGraph_TRUE))
    {
        // the structure of A is known to be symmetric
        Degree = G->rowdegree ;
    }
    else
    {
        // A is not known to be symmetric
        Degree = (byrow) ? G->rowdegree : G->coldegree ;
    }

    LG_ASSERT_MSG (Degree != NULL,
        LAGRAPH_PROPERTY_MISSING, "degree property unknown") ;

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    LG_TRY (LAGraph_Malloc ((void **) &samples, nsamples, sizeof (int64_t), msg)) ;

    //--------------------------------------------------------------------------
    // pick nsamples nodes at random and determine their degree
    //--------------------------------------------------------------------------

    // See also the hashed sampling method in LG_CC_FastSV6, which computes a
    // fast estimate of the mode of an integer vector.  This method does not
    // require a hash table.  However, the mode estimator in LG_CC_FastSV6
    // would be a good candidate to add as an LAGraph_SampleMode utility
    // function.

    GrB_Index n ;
    GRB_TRY (GrB_Vector_size (&n, Degree)) ;

    int64_t dsum = 0 ;
    for (int k = 0 ; k < nsamples ; k++)
    {
        uint64_t result = LG_Random60 (&seed) ;
        int64_t i = result % n ;
        // d = Degree (i)
        int64_t d ;
        GRB_TRY (GrB_Vector_extractElement (&d, Degree, i)) ;
        samples [k] = d ;
        dsum += d ;
    }

    // find the mean degree
    (*sample_mean) = ((double) dsum) / nsamples ;

    // find the median degree
    LG_qsort_1a (samples, nsamples) ;
    (*sample_median) = (double) samples [nsamples/2] ;

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_ALL ;
    return (GrB_SUCCESS) ;
}
