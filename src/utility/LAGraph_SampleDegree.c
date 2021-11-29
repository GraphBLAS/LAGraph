//------------------------------------------------------------------------------
// LAGraph_SampleDegree: sample the degree median and mean
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_SampleDegree computes estimates of the mean and median of the
// row or column degree of a graph.

#define LAGraph_FREE_ALL LAGraph_Free ((void **) &samples) ;

#include "LG_internal.h"

int LAGraph_SampleDegree        // returns 0 if successful, -1 if failure
(
    double *sample_mean,        // sampled mean degree
    double *sample_median,      // sampled median degree
    // input
    LAGraph_Graph G,        // graph of n nodes
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
    LG_CHECK (sample_mean == NULL, -1, "sample_mean is null") ;
    LG_CHECK (sample_median == NULL, -1, "sample_median is null") ;
    nsamples = LAGraph_MAX (nsamples, 1) ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;

    GrB_Vector Degree ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_structure_is_symmetric == LAGRAPH_TRUE))
    {
        // the structure of A is known to be symmetric
        Degree = G->rowdegree ;
    }
    else
    {
        // A is not known to be symmetric
        Degree = (byrow) ? G->rowdegree : G->coldegree ;
    }

    LG_CHECK (Degree == NULL, -1, "degree property unknown") ;

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    samples = LAGraph_Malloc (nsamples, sizeof (int64_t)) ;
    LG_CHECK (samples == NULL, -1, "out of memory") ;

    //--------------------------------------------------------------------------
    // pick nsamples nodes at random and determine their degree
    //--------------------------------------------------------------------------

    // See also the hashed sampling method in LG_CC_FastSV6, which computes a
    // fast estimate of the mode of an integer vector.  This method does not
    // require a hash table.  However, the mode estimator in LG_CC_FastSV6
    // would be a good candidate to add as an LAGraph_SampleMode utility
    // function.

    GrB_Index n ;
    GrB_TRY (GrB_Vector_size (&n, Degree)) ;

    int64_t dsum = 0 ;
    for (int k = 0 ; k < nsamples ; k++)
    {
        uint64_t result = LAGraph_Random60 (&seed) ;
        int64_t i = result % n ;
        // d = Degree (i)
        int64_t d ;
        GrB_TRY (GrB_Vector_extractElement (&d, Degree, i)) ;
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

    LAGraph_FREE_ALL ;
    return (0) ;
}
