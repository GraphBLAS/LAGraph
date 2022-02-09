//------------------------------------------------------------------------------
// LAGraph_DisplayGraph: print the contents of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DisplayGraph
(
    // input:
    LAGraph_Graph G,        // graph to display
    // TODO: use an enum for pr
    int pr,                 // 0: nothing, 1: terse, 2: summary, 3: all,
                            // 4: same as 2 but with %0.15g for doubles
                            // 5: same as 3 but with %0.15g for doubles
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
    pr = LAGraph_MAX (pr, -1) ;
    pr = LAGraph_MIN (pr, 5) ;

    //--------------------------------------------------------------------------
    // display the primary graph components
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;
    LAGraph_Kind kind = G->kind ;

    GrB_Index n, nvals ;
    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;
    char typename [LAGRAPH_MAX_NAME_LEN] ;
    char kindname [LAGRAPH_MAX_NAME_LEN] ;
    LG_TRY (LAGraph_MatrixTypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_KindName (kindname, kind, msg)) ;

    if (pr >= 0)
    {
        // print the basic scalar properties
        FPRINTF (f, "Graph: kind: %s, nodes: %g entries: %g type: %s\n",
            kindname, (double)n, (double)nvals, typename) ;
    }

    if (pr <= 0)
    {
        // nothing more to do
        return (GrB_SUCCESS) ;
    }

    // print the scalar cached properties
    FPRINTF (f, "  structural symmetry: ") ;
    switch (G->A_structure_is_symmetric)
    {
        case LAGRAPH_FALSE : FPRINTF (f, "unsymmetric") ; break ;
        case LAGRAPH_TRUE  : FPRINTF (f, "symmetric")   ; break ;
        default            : FPRINTF (f, "unknown")     ; break ;
    }
    if (G->ndiag >= 0) FPRINTF (f, "  self-edges: %g", (double) G->ndiag) ;
    FPRINTF (f, "\n") ;

    // pr = LAGraph_MAX (pr, 0) ;
    FPRINTF (f, "  adjacency matrix: ") ;

    LG_TRY (LAGraph_Matrix_print (A, pr, stdout, msg)) ;

    //--------------------------------------------------------------------------
    // display the cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix AT = G->AT ;
    if (AT != NULL)
    {
        FPRINTF (f, "  adjacency matrix transposed: ") ;
        LG_TRY (LAGraph_Matrix_print (AT, pr, stdout, msg)) ;
    }

    GrB_Vector rowdegree = G->rowdegree ;
    if (rowdegree != NULL)
    {
        FPRINTF (f, "  row degree: ") ;
        LG_TRY (LAGraph_Vector_print (rowdegree, pr, stdout, msg)) ;
    }

    GrB_Vector coldegree = G->coldegree ;
    if (coldegree != NULL)
    {
        FPRINTF (f, "  column degree: ") ;
        LG_TRY (LAGraph_Vector_print (coldegree, pr, stdout, msg)) ;
    }

    return (GrB_SUCCESS) ;
}
