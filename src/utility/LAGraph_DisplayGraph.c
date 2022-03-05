//------------------------------------------------------------------------------
// LAGraph_DisplayGraph: print the contents of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DisplayGraph
(
    // input:
    const LAGraph_Graph G,  // graph to display
    LAGraph_Print_Level pr, // print level (0 to 5)
    FILE *f,                // file to write to, must already be open
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear the msg and check the graph
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (f != NULL, GrB_NULL_POINTER) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    int prl = (int) pr ;
    prl = LAGRAPH_MAX (prl, 0) ;
    prl = LAGRAPH_MIN (prl, 5) ;
    if (prl == 0) return (GrB_SUCCESS) ;

    //--------------------------------------------------------------------------
    // display the primary graph components
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    LAGraph_Kind kind = G->kind ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    char typename [LAGRAPH_MAX_NAME_LEN] ;
    char kindname [LAGRAPH_MAX_NAME_LEN] ;
    LG_TRY (LAGraph_Matrix_TypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_KindName (kindname, kind, msg)) ;

    // print the basic scalar properties
    FPRINTF (f, "Graph: kind: %s, nodes: %g entries: %g type: %s\n",
        kindname, (double)n, (double)nvals, typename) ;

    // print the scalar cached properties
    FPRINTF (f, "  structural symmetry: ") ;
    switch (G->structure_is_symmetric)
    {
        case LAGraph_FALSE : FPRINTF (f, "unsymmetric") ; break ;
        case LAGraph_TRUE  : FPRINTF (f, "symmetric")   ; break ;
        default            : FPRINTF (f, "unknown")     ; break ;
    }
    if (G->ndiag >= 0) FPRINTF (f, "  self-edges: %g", (double) G->ndiag) ;
    FPRINTF (f, "\n") ;

    FPRINTF (f, "  adjacency matrix: ") ;

    LAGraph_Print_Level pr2 = (LAGraph_Print_Level) prl ;
    LG_TRY (LAGraph_Matrix_Print (A, pr2, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // display the cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix AT = G->AT ;
    if (AT != NULL)
    {
        FPRINTF (f, "  adjacency matrix transposed: ") ;
        LG_TRY (LAGraph_Matrix_Print (AT, pr2, stdout, msg)) ;
    }

    GrB_Vector rowdegree = G->rowdegree ;
    if (rowdegree != NULL)
    {
        FPRINTF (f, "  row degree: ") ;
        LG_TRY (LAGraph_Vector_Print (rowdegree, pr2, stdout, msg)) ;
    }

    GrB_Vector coldegree = G->coldegree ;
    if (coldegree != NULL)
    {
        FPRINTF (f, "  column degree: ") ;
        LG_TRY (LAGraph_Vector_Print (coldegree, pr2, stdout, msg)) ;
    }

    return (GrB_SUCCESS) ;
}
