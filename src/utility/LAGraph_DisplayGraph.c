//------------------------------------------------------------------------------
// LAGraph_DisplayGraph: print the contents of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DisplayGraph    // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to display
    int pr,                 // -1: nothing, 0: single line, 1: terse,
                            // 2: summary, 3: all,
                            // 4: same as 2 but with %0.15g for doubles
                            // 5: same as 3 but with %0.15g for doubles
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear the msg and check the graph
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LAGraph_TRY (LAGraph_CheckGraph (G, msg)) ;
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
    char *typename, *kindname ;
    LAGraph_TRY (LAGraph_TypeName (&typename, G->A_type, msg)) ;
    LAGraph_TRY (LAGraph_KindName (&kindname, kind, msg)) ;

    if (pr >= 0)
    {
        // print the basic scalar properties
        printf ("Graph: kind: %s, nodes: %ld entries: %ld type: %s\n",
            kindname, n, nvals, typename) ;

        // print the scalar cached properties
        printf ("    pattern symmetry: ") ;
        switch (G->A_pattern_is_symmetric)
        {
            case LAGRAPH_FALSE : printf ("unsymmetric") ; break ;
            case LAGRAPH_TRUE  : printf ("symmetric")   ; break ;
            default            : printf ("unknown")     ; break ;
        }
        if (G->ndiag >= 0) printf ("  self-edges: %ld", G->ndiag) ;
        printf ("\n") ;
    }

    pr = LAGraph_MAX (pr, 0) ;

    GrB_TRY (GxB_print (A, pr)) ;       // FIXME

    //--------------------------------------------------------------------------
    // display the cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix AT = G->AT ;
    if (AT != NULL)
    {
        GrB_TRY (GxB_print (AT, pr)) ;      // FIXME
    }

    GrB_Vector rowdegree = G->rowdegree ;
    if (rowdegree != NULL)
    {
        GrB_TRY (GxB_print (rowdegree, pr)) ;       // FIXME
    }

    GrB_Vector coldegree = G->coldegree ;
    if (coldegree != NULL)
    {
        GrB_TRY (GxB_print (coldegree, pr)) ;       // FIXME
    }

    return (0) ;
}
