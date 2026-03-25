//------------------------------------------------------------------------------
// LAGr_KatzCentrality: Katz centrality algorithm
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Karan Bhalla and Timothy A. Davis, Texas A&M University;

//------------------------------------------------------------------------------

#define LG_FREE_WORK                                                           \
{                                                                              \
    GrB_free(&x_prev);                                                         \
    GrB_free(&b);                                                              \
    GrB_free(&t);                                                              \
}

#define LG_FREE_ALL                                                            \
{                                                                              \
    LG_FREE_WORK;                                                              \
    GrB_free(&x);                                                              \
}

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_KatzCentrality
(
    // output:
    GrB_Vector *centrality,  
    int64_t *iters,   
    // input:
    LAGraph_Graph G,            
    double alpha,           
    double beta,            
    int64_t max_iter,           
    double tol,
    bool normalize,
    bool use_weights,               // uses edge weights if true, otherwise treats all edges as weight 1
    char* msg              
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;
    GrB_Index n = 0;
    GrB_Vector x = NULL, x_prev = NULL, b = NULL, t = NULL;

    LG_ASSERT (centrality != NULL && iters != NULL, GrB_NULL_POINTER) ;
    (*centrality) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    GrB_Matrix AT;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        AT = G->A ;
    }
    else
    {
        AT = G->AT ;
        LG_ASSERT_MSG (AT != NULL, LAGRAPH_NOT_CACHED, "G->AT is required") ;
    }

    // compute correct semiring based on whether edge weights are used
    GrB_Semiring semiring = use_weights ? GrB_PLUS_TIMES_SEMIRING_FP64 : GxB_PLUS_SECOND_FP64;

    if (use_weights)
    {
        // ensure min edge weight is available and nonnegative
        LG_TRY (LAGraph_Cached_EMin (G, msg)) ;
        LG_ASSERT_MSG (G->emin != NULL &&
            (G->emin_state == LAGraph_VALUE ||
             G->emin_state == LAGraph_BOUND),
            LAGRAPH_NOT_CACHED, "G->emin is required") ;

        double emin = 0 ;
        GRB_TRY (GrB_Scalar_extractElement_FP64 (&emin, G->emin)) ;
        LG_ASSERT_MSG (emin >= 0, GrB_INVALID_VALUE,
            "use_weights=true requires nonnegative edge weights") ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------
    
    GRB_TRY (GrB_Matrix_nrows (&n, AT)) ;

    GRB_TRY(GrB_Vector_new(&x, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&x_prev, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&b, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&t, GrB_FP64, n));

    GRB_TRY(GrB_assign(x, NULL, NULL, 0.0, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(x_prev, NULL, NULL, 0.0, GrB_ALL, n, NULL));
    GRB_TRY(GrB_assign(b, NULL, NULL, beta, GrB_ALL, n, NULL));

    // first iteration is always done
    double rdiff = 1 ;   

    // TODO: determine best way to stop 
    for ((*iters) = 0 ; ; (*iters)++)
    {
        // check for convergence failure
        LG_ASSERT_MSGF ((*iters) < max_iter, LAGRAPH_CONVERGENCE_FAILURE,
            "katz centrality failed to converge in %d iterations", max_iter) ;

        // swap x and x_prev
        GrB_Vector temp = x_prev ; x_prev = x ; x = temp ;

        // x = A' * x_prev
        GRB_TRY (GrB_mxv (x, NULL, NULL, semiring, AT, x_prev, NULL)) ;

        // x = alpha * x + beta
        GRB_TRY (GrB_apply (x, NULL, NULL, GrB_TIMES_FP64, alpha, x, NULL)) ;
        GRB_TRY (GrB_eWiseAdd (x, NULL, NULL, GrB_PLUS_FP64, x, b, NULL)) ;

        // t = x - x_prev
        GRB_TRY (GrB_eWiseAdd (t, NULL, NULL, GrB_MINUS_FP64, x, x_prev, NULL)) ;
        // t = abs (t)
        GRB_TRY (GrB_apply (t, NULL, NULL, GrB_ABS_FP64, t, NULL)) ;
        // rdiff = sum (t)
        GRB_TRY (GrB_reduce (&rdiff, NULL, GrB_PLUS_MONOID_FP64, t, NULL)) ;

        if (rdiff < tol) break ;        // TODO: revisit
    }

    // normalize using the L2 norm if flag is set
    if (normalize)
    {
        double sumsq = 0 ;

        // sumsq = x' * x = sum (x .* x)
        GRB_TRY (GrB_eWiseMult (t, NULL, NULL, GrB_TIMES_FP64, x, x, NULL)) ;
        GRB_TRY (GrB_reduce (&sumsq, NULL, GrB_PLUS_MONOID_FP64, t, NULL)) ;

        if (sumsq > 0)
        {
            GRB_TRY (GrB_apply (x, NULL, NULL, GrB_DIV_FP64, x, sqrt(sumsq), NULL)) ;
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*centrality) = x;
    LG_FREE_WORK;

    return GrB_SUCCESS;
}
