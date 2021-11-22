//------------------------------------------------------------------------------
// LAGraph_cc_boruvka.c
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

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_FREE_WORK ;             \
    GrB_free (&f) ;                 \
}

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

#include "LAGraphX.h"
#include "LG_internal.h"

//****************************************************************************
// w[index[i]] = min(w[index[i]], s[i]) for i in [0..n-1]
static GrB_Info Reduce_assign (GrB_Vector w,
                               GrB_Vector s,
                               GrB_Index *index,
                               GrB_Index n)
{
    GrB_Index *mem = (GrB_Index *) LAGraph_Malloc (3*n, sizeof (GrB_Index)) ;
    // FIXME: check for malloc failure
    GrB_Index *ind = mem, *sval = mem + n, *wval = sval + n;
    GrB_Vector_extractTuples(ind, wval, &n, w) ;
    GrB_Vector_extractTuples(ind, sval, &n, s) ;
    for (GrB_Index i = 0; i < n; i++)
    {
        if (sval[i] < wval[index[i]])
        {
            wval[index[i]] = sval[i];
        }
    }
    GrB_Vector_clear(w) ;
    GrB_Vector_build(w, ind, wval, n, GrB_PLUS_UINT64) ;
    LAGraph_Free ((void **) &mem) ;
    return GrB_SUCCESS;
}

//------------------------------------------------------------------------------
// select_func: IndexUnaryOp for pruning entries from S
//------------------------------------------------------------------------------

// FIXME: this uses global variables; fix this 

static GrB_Index *I, *V ;

#if 1
void select_func (void *z, const void *x, 
                 const GrB_Index i, const GrB_Index j,
                 const void *y)
{
    (*(bool *) z) = (V [i] != V [j]) ;
}
#else
bool select_func (const GrB_Index i, const GrB_Index j,
                 const void *x, const void *thunk)
{
    return (V [i] != V [j]) ;
}
#endif

//****************************************************************************

GrB_Info LAGraph_cc_boruvka
(
    GrB_Vector *result,     // output: array of component identifiers
    GrB_Matrix A,           // input matrix
    bool sanitize           // if true, ensure A is symmetric
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    printf ("------------------------ Boruvka\n") ;
    GrB_Info info ;
    char *msg = NULL ;
    if (result == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n;
    GrB_Vector f = NULL, gp = NULL, mnp = NULL, ccmn = NULL, ramp = NULL,
        mask = NULL ;
    #if 1
    GrB_IndexUnaryOp select_op = NULL ;
    #else
    GxB_SelectOp select_op = NULL ;
    #endif

    GrB_TRY (GrB_Matrix_nrows (&n, A)) ;

    // FIXME: make S boolean and iso-valued
    GrB_Matrix S ;
    if (sanitize)
    {
        GrB_Descriptor desc;
        GrB_TRY (GrB_Matrix_new (&S, GrB_BOOL, n, n)) ;
        GrB_TRY (GrB_eWiseAdd (S, NULL, NULL, GrB_LOR, A, A, GrB_DESC_T1)) ;
    }
    else
    {
        // Use the input as-is, and assume it is binary and symmetric
        GrB_TRY (GrB_Matrix_dup (&S, A)) ;
    }

    GrB_TRY (GrB_Vector_new (&f, GrB_UINT64, n)) ;      //
    GrB_TRY (GrB_Vector_new (&gp, GrB_UINT64, n)) ;     // grandparents
    GrB_TRY (GrB_Vector_new (&mnp, GrB_UINT64, n)) ;
    GrB_TRY (GrB_Vector_new (&ccmn, GrB_UINT64, n)) ;
    GrB_TRY (GrB_Vector_new (&ramp, GrB_UINT64, n)) ;
    GrB_TRY (GrB_Vector_new (&mask, GrB_BOOL, n)) ;

    // prepare
    I = (GrB_Index*) LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    V = (GrB_Index*) LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    // FIXME: use ROWINDEX operator
    for (GrB_Index i = 0; i < n; i++)
    {
        I[i] = V[i] = i;
    }
    GrB_TRY (GrB_Vector_build (f, I, V, n, GrB_PLUS_UINT64)) ;
    GrB_TRY (GrB_assign (ramp, NULL, NULL, f, GrB_ALL, n, NULL)) ;

    #if 1
    GrB_TRY (GrB_IndexUnaryOp_new (&select_op, select_func,
        GrB_BOOL, /* aij: ignored */ GrB_BOOL, /* y: ignored */ GrB_BOOL)) ;
    #else
    GrB_TRY (GxB_SelectOp_new (&select_op, select_func, NULL, NULL)) ;
    #endif

    GrB_Index nvals ;
    GrB_TRY (GrB_Matrix_nvals (&nvals, S)) ;
    double tselect = 0 ;
    printf ("S nvals %ld\n", nvals) ;

    //--------------------------------------------------------------------------
    // find the connected components
    //--------------------------------------------------------------------------

    int64_t niters ;
    for (niters = 0 ; nvals > 0 ; niters++)
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

        // gb = f (f)
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
            // gb = f (f)
            GrB_TRY (GrB_Vector_extractTuples (I, V, &n, f)) ;
            GrB_TRY (GrB_extract (gp, NULL, NULL, f, V, n, NULL)) ;

            // mask = (f != gp)
            GrB_TRY (GrB_eWiseMult (mask, NULL, NULL, GrB_NE_UINT64, f, gp,
                NULL)) ;

            // swap f and gp
            GrB_Vector t = f ; f = gp ; gp = t ;

            // diff = or (mask)
            GrB_TRY (GrB_reduce (&diff, NULL, GrB_LOR_MONOID_BOOL, mask, NULL)) ;
        }

        //----------------------------------------------------------------------
        // remove the edges inside each connected component
        //----------------------------------------------------------------------

        double ttt = omp_get_wtime ( ) ;
        #if 1
        GrB_TRY (GrB_select (S, NULL, NULL, select_op, S, (bool) false, NULL)) ;
        #else
        GrB_TRY (GxB_select (S, NULL, NULL, select_op, S, 0, NULL)) ;
        #endif
        GrB_TRY (GrB_Matrix_nvals (&nvals, S)) ;
        ttt = omp_get_wtime ( ) - ttt ;
        printf ("S nvals %ld select: %g\n", nvals, ttt) ;
        tselect += ttt ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    printf ("Boruvka iters %ld select time %g\n", niters, tselect) ;
    (*result) = f ;
    LAGraph_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
