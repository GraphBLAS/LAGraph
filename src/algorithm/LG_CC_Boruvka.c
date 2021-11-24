//------------------------------------------------------------------------------
// LG_CC_Boruvka.c
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// Code is based on Boruvka's minimum spanning forest algorithm.
// Contributed by Yongzhe Zhang (zyz915@gmail.com).
// revised by Tim Davis (davis@tamu.edu)

#include "LG_internal.h"

//------------------------------------------------------------------------------
// Reduce_assign
//------------------------------------------------------------------------------

// w[index[i]] = min(w[index[i]], s[i]) for i in [0..n-1]

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_Free ((void **) &mem) ; \
}

static GrB_Info Reduce_assign
(
    GrB_Vector w,       // input/output vector of size n
    GrB_Vector s,       // input vector of size n
    GrB_Index *index,   // index array of size n
    GrB_Index n
)
{
    char *msg = NULL ;
    GrB_Index *mem = (GrB_Index *) LAGraph_Malloc (3*n, sizeof (GrB_Index)) ;
    LG_CHECK (mem == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;
    GrB_Index *ind = mem, *sval = mem + n, *wval = sval + n ;
    GrB_TRY (GrB_Vector_extractTuples (ind, wval, &n, w)) ;
    GrB_TRY (GrB_Vector_extractTuples (ind, sval, &n, s)) ;
    for (GrB_Index i = 0 ; i < n ; i++)
    {
        if (sval [i] < wval [index [i]])
        {
            wval [index [i]] = sval [i] ;
        }
    }
    GrB_TRY (GrB_Vector_clear (w)) ;
    GrB_TRY (GrB_Vector_build (w, ind, wval, n, GrB_PLUS_UINT64)) ;
    LAGraph_FREE_ALL ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// select_func: IndexUnaryOp for pruning entries from S
//------------------------------------------------------------------------------

// TODO: this uses global variables; fix this.  Use y to pass in 2 pointers.

static GrB_Index *V ;

void select_func (void *z, const void *x, 
                 const GrB_Index i, const GrB_Index j, const void *y)
{
    (*((bool *) z)) = (V [i] != V [j]) ;
}

//------------------------------------------------------------------------------
// LG_CC_Boruvka
//------------------------------------------------------------------------------

#undef  LAGraph_FREE_ALL
#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_FREE_WORK ;             \
    GrB_free (&f) ;                 \
}

#undef  LAGraph_FREE_WORK
#define LAGraph_FREE_WORK           \
{                                   \
    LAGraph_Free ((void **) &I) ;   \
    LAGraph_Free ((void **) &V) ;   \
    GrB_free (&gp) ;                \
    GrB_free (&mnp) ;               \
    GrB_free (&ccmn) ;              \
    GrB_free (&ramp) ;              \
    GrB_free (&mask) ;              \
    GrB_free (&select_op) ;         \
}

int LG_CC_Boruvka
(
    // output
    GrB_Vector *component,  // output: array of component identifiers
    // inputs
    LAGraph_Graph G,        // input graph, not modified
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Info info ;
    GrB_Index n, *I = NULL ;
    GrB_Vector f = NULL, gp = NULL, mnp = NULL, ccmn = NULL, ramp = NULL,
        mask = NULL ;
    GrB_IndexUnaryOp select_op = NULL ;
    GrB_Matrix S = NULL ;
    V = NULL ;

    LG_CLEAR_MSG ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;
    LG_CHECK (component == NULL, GrB_NULL_POINTER, "input is NULL") ;

    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_structure_is_symmetric == LAGRAPH_TRUE))
    {
        // A must be symmetric
        ;
    }
    else
    {
        // A must not be unsymmetric
        LG_CHECK (false, -1, "input must be symmetric") ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    // S = structure of G->A
    LAGraph_TRY (LAGraph_Structure (&S, G->A, msg)) ;

    GrB_TRY (GrB_Matrix_nrows (&n, S)) ;
    GrB_TRY (GrB_Vector_new (&f, GrB_UINT64, n)) ;      // final result
    GrB_TRY (GrB_Vector_new (&gp, GrB_UINT64, n)) ;     // grandparents
    GrB_TRY (GrB_Vector_new (&mnp, GrB_UINT64, n)) ;    // min neighbor parent
    GrB_TRY (GrB_Vector_new (&ccmn, GrB_UINT64, n)) ;   // cc's min neighbor
    GrB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n)) ;     // various uses

    V = (GrB_Index*) LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    LG_CHECK (V == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;
    #if !LG_SUITESPARSE
    // I is not needed for SuiteSparse and remains NULL
    I = (GrB_Index*) LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    LG_CHECK (I == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;
    #endif

    // f = 0:n-1, and copy to ramp
    GrB_TRY (GrB_assign (f, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GrB_TRY (GrB_apply  (f, NULL, NULL, GrB_ROWINDEX_INT64, f, 0, NULL)) ;
    GrB_TRY (GrB_Vector_dup (&ramp, f)) ;
    GrB_TRY (GrB_Vector_extractTuples (I, V, &n, f)) ;

    GrB_TRY (GrB_IndexUnaryOp_new (&select_op, select_func,
        GrB_BOOL, /* aij: ignored */ GrB_BOOL, /* y: ignored */ GrB_BOOL)) ;

    GrB_Index nvals ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, S)) ;

    //--------------------------------------------------------------------------
    // find the connected components
    //--------------------------------------------------------------------------

    while (nvals > 0)
    {

        //----------------------------------------------------------------------
        // mnp[u] = u's minimum neighbor's parent for all nodes u
        //----------------------------------------------------------------------

        // every vertex points to a root vertex at the begining
        GrB_TRY (GrB_assign (mnp, NULL, NULL, n, GrB_ALL, n, NULL)) ;
        GrB_TRY (GrB_mxv (mnp, NULL, GrB_MIN_UINT64,
                    GrB_MIN_SECOND_SEMIRING_UINT64, S, f, NULL)) ;

        //----------------------------------------------------------------------
        // find the minimum neighbor
        //----------------------------------------------------------------------

        // ccmn[u] = connect component's minimum neighbor | if u is a root
        //         = n                                    | otherwise
        GrB_TRY (GrB_assign (ccmn, NULL, NULL, n, GrB_ALL, n, NULL)) ;
        GrB_TRY (Reduce_assign (ccmn, mnp, V, n)) ;

        //----------------------------------------------------------------------
        // f[u] = ccmn[u] if ccmn[u] != n
        //----------------------------------------------------------------------

        // mask = (ccnm != n)
        GrB_TRY (GrB_apply (mask, NULL, NULL, GrB_NE_UINT64, ccmn, n, NULL)) ;
        // f<mask> = ccmn
        GrB_TRY (GrB_assign (f, mask, NULL, ccmn, GrB_ALL, n, NULL)) ;

        //----------------------------------------------------------------------
        // select new roots
        //----------------------------------------------------------------------

        // identify all the vertex pairs (u, v) where f[u] == v and f[v] == u
        // and then select the minimum of u, v as the new root;
        // if (f[f[i]] == i) f[i] = min(f[i], i)

        // compute grandparents: gp = f (f)
        GrB_TRY (GrB_Vector_extractTuples (I, V, &n, f)) ;
        GrB_TRY (GrB_extract (gp, NULL, NULL, f, V, n, NULL)) ;

        // mask = (gp == 0:n-1)
        GrB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_EQ_UINT64, gp, ramp,
            NULL)) ;
        // f<mask> = min (f, ramp)
        GrB_TRY (GrB_assign (f, mask, GrB_MIN_UINT64, ramp, GrB_ALL, n, NULL)) ;

        //----------------------------------------------------------------------
        // shortcutting f[i] = f[f[i]] until f does not change
        //----------------------------------------------------------------------

        bool diff = true ;
        while (diff)
        {
            // compute grandparents: gp = f (f)
            GrB_TRY (GrB_Vector_extractTuples (I, V, &n, f)) ;
            GrB_TRY (GrB_extract (gp, NULL, NULL, f, V, n, NULL)) ;

            // mask = (f != gp)
            GrB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_NE_UINT64, f, gp,
                NULL)) ;

            // swap f and gp
            GrB_Vector t = f ; f = gp ; gp = t ;

            // diff = or (mask)
            GrB_TRY (GrB_reduce (&diff, NULL, GrB_LOR_MONOID_BOOL, mask,
                NULL)) ;
        }

        //----------------------------------------------------------------------
        // remove the edges inside each connected component
        //----------------------------------------------------------------------

        GrB_TRY (GrB_select (S, NULL, NULL, select_op, S, (bool) false, NULL)) ;
        GrB_TRY (GrB_Matrix_nvals (&nvals, S)) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*component) = f ;
    LAGraph_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
