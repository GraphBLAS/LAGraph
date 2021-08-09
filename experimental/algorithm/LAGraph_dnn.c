//------------------------------------------------------------------------------
// LAGraph_dnn: sparse deep neural network
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// LAGraph_dnn: sparse deep neural network.  Contributed by Tim Davis,
// Texas A&M University.  Based on inferenceReLUvec.m by Jeremy Kepner, MIT.

// Performs ReLU inference using input feature vectors Y0.

// See http://graphchallenge.org/ for a description of the algorithm.

// On input, Y0 is the initial feature vectors, of size nfeatures-by-nneurons.
// This format uses the graph convention that A(i,j) is the edge (i,j),
// Each row of Y0 is a single feature.

// W is an array of size nlayers of sparse matrices.  Each W[layer] matrix has
// the same size: nneurons-by-nneurons.  W[layer] represents the DNN weights
// for that layer.

// The Bias[layer] matrices are diagaonal, and the same size as W[layer].

// All matrices must have the same type:  either GrB_FP32 or GrB_FP64.

// On output, Y is the computed result, of the same size and type as Y0.

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&gt0) ;       \
    GrB_free (&ymax) ;      \
    GrB_free (&M) ;         \
    GrB_free (Yhandle) ;    \
}

#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING
#include <LAGraph.h>
#include <LAGraphX.h>

//****************************************************************************
#define F_UNARY(f)  ((void (*)(void *, const void *)) f)

// unary ops to check if greater than zero
void LAGraph_gt0_fp32
(
    bool *z,
    const float *x
)
{
    (*z) = ((*x) > 0) ;
}

void LAGraph_gt0_fp64
(
    bool *z,
    const double *x
)
{
    (*z) = ((*x) > 0) ;
}

// unary operators to threshold a max value for DNN
void LAGraph_ymax_fp32
(
    float *z,
    const float *x
)
{
    (*z) = fminf ((*x), (float) 32.0) ;
}

void LAGraph_ymax_fp64
(
    double *z,
    const double *x
)
{
    (*z) = fmin ((*x), (double) 32.0) ;
}

//****************************************************************************
GrB_Info LAGraph_dnn    // returns GrB_SUCCESS if successful
(
    // output
    GrB_Matrix *Yhandle,    // Y, created on output
    // input: not modified
    GrB_Matrix *W,      // W [0..nlayers-1], each nneurons-by-nneurons
    GrB_Matrix *Bias,   // Bias [0..nlayers-1], diagonal nneurons-by-nneurons
    int nlayers,        // # of layers
    GrB_Matrix Y0       // input features: nfeatures-by-nneurons
)
{
    GrB_Info info ;

#if !defined(LG_SUITESPARSE)
    // GxB_type, semirings and select required
    return GrB_PANIC;
#else
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (Yhandle == NULL || W == NULL || Bias == NULL || Y0 == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    // unary ops to check if greater than zero
    GrB_Matrix Y = NULL ;
    GrB_Matrix M = NULL ;
    (*Yhandle) = NULL ;

    GrB_UnaryOp gt0 = NULL ;
    GrB_UnaryOp ymax = NULL ;

    //--------------------------------------------------------------------------
    // create the unary greater-than-zero operators
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // create the unary YMAX operators
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    GrB_Semiring plus_times, plus_plus ;
    GrB_UnaryOp id ;

    GrB_Type type ;
    LAGRAPH_OK (GxB_Matrix_type (&type, Y0)) ;

    if (type == GrB_FP32)
    {
        LAGRAPH_OK (GrB_UnaryOp_new (&gt0,
                                     F_UNARY (LAGraph_gt0_fp32),
                                     GrB_BOOL, GrB_FP32)) ;

        LAGRAPH_OK (GrB_UnaryOp_new (&ymax,
                                     F_UNARY (LAGraph_ymax_fp32),
                                     GrB_FP32, GrB_FP32)) ;

        plus_times = GrB_PLUS_TIMES_SEMIRING_FP32 ;
        plus_plus  = GxB_PLUS_PLUS_FP32 ;
        id         = GrB_IDENTITY_FP32 ;
    }
    else if (type == GrB_FP64)
    {
        LAGRAPH_OK (GrB_UnaryOp_new (&gt0,
                                     F_UNARY (LAGraph_gt0_fp64),
                                     GrB_BOOL, GrB_FP64)) ;

        LAGRAPH_OK (GrB_UnaryOp_new (&ymax,
                                     F_UNARY (LAGraph_ymax_fp64),
                                     GrB_FP64, GrB_FP64)) ;

        plus_times = GrB_PLUS_TIMES_SEMIRING_FP64 ;
        plus_plus  = GxB_PLUS_PLUS_FP64 ;
        id         = GrB_IDENTITY_FP64 ;
    }
    else
    {
        return (GrB_DOMAIN_MISMATCH) ;
    }

    for (int layer = 0 ; layer < nlayers ; layer++)
    {
        GrB_Type type2 ;
        LAGRAPH_OK (GxB_Matrix_type (&type2, W [layer])) ;
        if (type != type2) return (GrB_DOMAIN_MISMATCH) ;
        LAGRAPH_OK (GxB_Matrix_type (&type2, Bias [layer])) ;
        if (type != type2) return (GrB_DOMAIN_MISMATCH) ;
    }

    //--------------------------------------------------------------------------
    // create the output matrix Y
    //--------------------------------------------------------------------------

    GrB_Index nfeatures, nneurons ;
    LAGRAPH_OK (GrB_Matrix_nrows (&nfeatures, Y0)) ;
    LAGRAPH_OK (GrB_Matrix_ncols (&nneurons,  Y0)) ;
    LAGRAPH_OK (GrB_Matrix_new (&Y, type, nfeatures, nneurons)) ;
    LAGRAPH_OK (GrB_Matrix_new (&M, GrB_BOOL, nfeatures, nneurons)) ;

    //--------------------------------------------------------------------------
    // propagate the features through the neuron layers
    //--------------------------------------------------------------------------

    for (int layer = 0 ; layer < nlayers ; layer++)
    {
        // Y = Y * W [layer], using the conventional PLUS_TIMES semiring
        LAGRAPH_OK (GrB_mxm (Y, NULL, NULL, plus_times,
            ((layer == 0) ? Y0 : Y), W [layer], NULL)) ;

        // Y = Y * Bias [layer], using the PLUS_PLUS semiring.  This computes
        // Y(i,j) += Bias [layer] (j,j) for each entry Y(i,j).  It does not
        // introduce any new entries in Y.
        LAGRAPH_OK (GrB_mxm (Y, NULL, NULL, plus_plus, Y, Bias [layer], NULL)) ;

        // delete entries from Y: keep only those entries greater than zero
        #if LG_SUITESPARSE
        // using SuiteSparse:GraphBLAS 3.0.0 or later.
        LAGRAPH_OK (GxB_select (Y, NULL, NULL, GxB_GT_ZERO, Y, NULL, NULL)) ;

        #else
        // using SuiteSparse v2.x or earlier, or any other GraphBLAS library.
        LAGRAPH_OK (GrB_apply (M, NULL, NULL, gt0, Y, NULL)) ;
        LAGRAPH_OK (GrB_apply (Y, M, NULL, id, Y, GrB_DESC_R)) ;
        #endif

        // threshold maximum values: Y (Y > 32) = 32
        // t = omp_get_wtime ( ) ;
        LAGRAPH_OK (GrB_apply (Y, NULL, NULL, ymax, Y, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    GrB_free (&M) ;
    (*Yhandle) = Y ;
    return (GrB_SUCCESS) ;
#endif
}
