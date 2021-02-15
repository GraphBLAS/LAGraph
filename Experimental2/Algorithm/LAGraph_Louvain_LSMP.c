//------------------------------------------------------------------------------
// LAGraph_Louvain_LSMP: Louvain Community Detection
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Based on the paper Linear Algebraic Louvain Method in Python
// by Tze Meng Low, Daniele G. Spampinato, Scott McMillan and Michel Pelletier (LSMP)
// Translated to from Python to C by Michel Pelletier
// Python Notebook: https://github.com/Graphegon/pygraphblas/blob/main/demo/Louvain.ipynb

//------------------------------------------------------------------------------

#include "LG_internal.h"

#define LAGRAPH_FREE_WORK          \
{                                  \
    GrB_free (&empty) ;             \
    GrB_free (&Sj) ;               \
    GrB_free (&v) ;                \
    GrB_free (&q) ;                \
    GrB_free (&t) ;                \
    GrB_free (&ApAT) ;             \
    GrB_free (&S) ;                \
}

#define LAGRAPH_FREE_ALL           \
{                                  \
    LAGRAPH_FREE_WORK ;            \
    GrB_free (&result) ;           \
}

int LAGraph_Louvain_LSMP // returns -1 on failure, 0 on success
(
    // outputs:
    GrB_Vector *community_assignment, // community assignment for node i
    // inputs:
    LAGraph_Graph G,        // input graph
    int itermax,            // maximum number of iterations (typically 100)
    int *iters,             // output: number of iterations taken
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector result = NULL, empty = NULL, Sj = NULL, k = NULL, v = NULL, q = NULL, t = NULL ;
    GrB_Matrix ApAT = NULL, S = NULL, W = NULL, GA = NULL;
    GxB_Scalar max;

    LG_CHECK (community_assignment == NULL, -1, "community_assignment is NULL") ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    GrB_Index n, kn, tn, r;
    float *tx;
    (*community_assignment) = NULL ;
    bool changed = true;
    float m, f, max_f;
    GrB_Index *ks, *ts;
    GrB_Index kj;

    GA = G->A ;

    GrB_TRY (GrB_Matrix_nrows (&n, GA)) ;

    GrB_TRY (GrB_Matrix_new (&S,     GrB_BOOL, n, n)) ;
    GrB_TRY (GrB_Matrix_new (&ApAT,  GrB_FP32, n, n)) ;
    GrB_TRY (GrB_Vector_new (&Sj,    GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&k,     GrB_FP32, n)) ;
    GrB_TRY (GrB_Vector_new (&empty, GrB_BOOL, n)) ;
    GrB_TRY (GrB_Vector_new (&v,     GrB_FP32, n)) ;
    GrB_TRY (GrB_Vector_new (&q,     GrB_FP32, n)) ;
    GrB_TRY (GrB_Vector_new (&t,     GrB_FP32, n)) ;
    GrB_TRY (GxB_Scalar_new (&max,   GrB_FP32)) ;

    // ApAT = A.T + A
    GrB_TRY (GrB_eWiseAdd  (ApAT, NULL, NULL, GrB_PLUS_FP32, GA, GA, GrB_DESC_T0)) ;

    GrB_TRY (GrB_reduce (k,    NULL, NULL, GrB_PLUS_MONOID_FP32, GA, NULL)) ; // k = A.reduce_vector()
    GrB_TRY (GrB_reduce (&m,   NULL, GrB_PLUS_MONOID_FP32, k, NULL)) ;       // m = k.reduce_int()
    GrB_TRY (GrB_apply  (k,    NULL, NULL, GrB_AINV_FP32, k, NULL)) ;              // k = (-k) / m
    GrB_TRY (GrB_apply  (k,    NULL, NULL, GrB_DIV_FP32, k, m, NULL)) ;

    GrB_TRY (GrB_Vector_size (&kn, k)) ;
    ks = LAGraph_Malloc (kn, sizeof (GrB_Index)) ;
    GrB_TRY (GrB_Vector_extractTuples (ks, (double*)NULL, &kn, k)) ;

    //     while changed and i < max_iters:
    for ((*iters) = 0 ; (*iters) < itermax && changed ; (*iters)++)
    {
        changed = false ;

        // for j in set(k.indexes):
        for (GrB_Index j = 0; j < kn; j++) {
            kj = ks[j];
            // Sj = S[j]
            GrB_TRY (GrB_extract  (Sj,     NULL, NULL, S, GrB_ALL, n, kj, GrB_DESC_T0)) ;
            // S[j] = empty
            GrB_TRY (GrB_assign   (S,      NULL, NULL, empty, kj, GrB_ALL, n, NULL)) ;
            GrB_TRY (GrB_extract  (v,      NULL, NULL, GA, GrB_ALL, n, kj, GrB_DESC_T0)) ;
            // v += k[j]
            GrB_TRY (GrB_Vector_extractElement    (&f,     v, kj)) ;
            GrB_TRY (GrB_apply    (v,      NULL, NULL, GrB_PLUS_FP32, v, f, NULL)) ;
            // q = v @ S
            GrB_TRY (GrB_vxm      (q,      NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP32, v, S, NULL)) ;
            // t = q.select('max')
            GrB_TRY (GrB_reduce   (&max_f, NULL, GrB_MAX_MONOID_FP32, q, NULL)) ;
            GrB_TRY (GxB_Scalar_setElement    (max,    max_f)) ;
            GrB_TRY (GxB_select   (t,      NULL, NULL, GxB_EQ_THUNK, q, max, NULL)) ;

            // if t:
            GrB_TRY (GrB_Vector_nvals (&tn, t)) ;
            if (tn == 0)
                continue;
             
            ts = LAGraph_Malloc (tn, sizeof (GrB_Index));
            GrB_TRY (GrB_Vector_extractTuples (ts, (float*)NULL, &tn, t)) ;

            if (tn == 1)
                r = ts[0] ;
            else
                r = ts[rand() % tn] ;

            GrB_TRY (GrB_Matrix_setElement(S, true, j, r)) ; //  S[j, r] = True
            if (GrB_Vector_extractElement(&changed, Sj, r)
                == GrB_NO_VALUE)                             //  if Sj.get(r) is None:
                changed = true;
			LAGraph_FREE (ts) ;
        }
    }

    // S.cast(INT64).apply(INT64.POSITIONJ).reduce_vector()
    GrB_TRY (GrB_Matrix_new (&W,      GrB_INT64, n, n)) ;
    GrB_TRY (GrB_apply      (W,       NULL, NULL, GxB_POSITIONJ_INT64, S, NULL)) ;
    GrB_TRY (GrB_Vector_new (&result, GrB_INT64, n)) ;
    GrB_TRY (GrB_reduce     (result,  NULL, NULL, GxB_ANY_INT64, W, NULL)) ;

    (*community_assignment) = result ;
    LAGRAPH_FREE_WORK ;
    return (0) ;
}
