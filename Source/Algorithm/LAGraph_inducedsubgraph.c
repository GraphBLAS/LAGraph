//------------------------------------------------------------------------------
// LAGraph_inducedsubgraph: extract the induced subgraph for a set of vertices
//------------------------------------------------------------------------------

/*
    LAGraph:  graph algorithms based on GraphBLAS

    Copyright 2019 LAGraph Contributors.

    (see Contributors.txt for a full list of Contributors; see
    ContributionInstructions.txt for information on how you can Contribute to
    this project).

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// LAGraph_inducedsubgraph: Contributed by Gabor Szarnyas.
// Budapest University of Technology and Economics
// (with accented characters: G\'{a}bor Sz\'{a}rnyas).

// Compute the

#include "LAGraph_internal.h"

#define LAGRAPH_FREE_ALL            \
{                                   \
    LAGRAPH_FREE (C) ;              \
    LAGRAPH_FREE (type) ;           \
}

typedef struct
{
    const GrB_Index nv; // number of vertices
    const GrB_Index* V; // vertex indices
} V_ind_struct;

int cmpfunc(const void *a, const void *b) {
    return ( *(GrB_Index *)a - *(GrB_Index *)b );
}

static bool func(const GrB_Index i, const GrB_Index j, const GrB_Index nrows,
                 const GrB_Index ncols, const void *x, const void *thunk)
{
    V_ind_struct* indices = (V_ind_struct*) (thunk);

    if (bsearch (&i, indices->V, indices->nv, sizeof (GrB_Index), cmpfunc))
    {
        return false;
    }
    if (bsearch (&j, indices->V, indices->nv, sizeof (GrB_Index), cmpfunc))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------

GrB_Info LAGraph_inducedsubgraph // compute the subgraph induced by vertices V
                                 // in A
(
    GrB_Matrix *Chandle,         // output matrix
    const GrB_Matrix A,          // input matrix
    const GrB_Index *V,          // sorted list of vertex indices
    GrB_Index nv,                // number of vertex indices
    bool use_select              // use GxB_select
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;

    GrB_Type type;
    GrB_Index n ;
    GrB_Matrix C = NULL ;

    LAGRAPH_OK (GxB_Matrix_type (&type, A)) ;
    LAGRAPH_OK (GrB_Matrix_nrows (&n, A)) ;
    LAGRAPH_OK (GrB_Matrix_new (&C, type, n, n)) ;

    if (use_select)
    {
        V_ind_struct indices = {.nv = nv, .V = V};

        GrB_Type V_ind_type;
        LAGRAPH_OK (GrB_Type_new (&V_ind_type, sizeof(indices)));

        GxB_Scalar Thunk;
        LAGRAPH_OK (GxB_Scalar_new (&Thunk, V_ind_type)) ;
        LAGRAPH_OK (GxB_Scalar_setElement_UDT (Thunk, &indices)) ;

        GxB_SelectOp sel_op;
        LAGRAPH_OK (GxB_SelectOp_new (&sel_op, func, GrB_NULL, V_ind_type));
        LAGRAPH_OK (GxB_select (C, NULL, NULL, sel_op, A, Thunk, NULL)) ;

        LAGRAPH_FREE(sel_op)
        LAGRAPH_FREE(Thunk)
        LAGRAPH_FREE(V_ind_type)
    }
    else
    {
        GrB_Matrix D; // diagonal matrix used to select rows/columns

        LAGRAPH_OK (GrB_Matrix_new(&D, GrB_BOOL, n, n));

        bool* X = LAGraph_malloc (nv, sizeof (GrB_BOOL)) ;
        if (X == NULL)
        {
            LAGRAPH_ERROR ("out of memory", GrB_OUT_OF_MEMORY) ;
        }
        for (GrB_Index i = 0; i < nv; i++) {
            X[i] = true;
        }
        LAGRAPH_OK (GrB_Matrix_build (D, V, V, X, nv, GrB_LOR)) ;

        LAGRAPH_OK (GrB_mxm (C, NULL, NULL, GxB_ANY_SECOND_FP64, D, A, NULL)) ;
        LAGRAPH_OK (GrB_mxm (C, NULL, NULL, GxB_ANY_FIRST_FP64,  C, D, NULL)) ;

        LAGRAPH_FREE(D);
    }
    (*Chandle) = C ;

    return (GrB_SUCCESS) ;
}

