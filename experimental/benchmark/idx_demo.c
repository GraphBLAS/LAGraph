//------------------------------------------------------------------------------
// LAGraph/src/benchmark/idx_demo.c: benchmark for user-defined IndexUnaryOp
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause

//------------------------------------------------------------------------------

// Usage:  tcc_demo < matrixmarketfile.mtx
//         tcc_demo matrixmarketfile.mtx
//         tcc_demo matrixmarketfile.grb

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"

#define LAGraph_FREE_ALL            \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    if (Ap != NULL) free (Ap) ;     \
    if (Ai != NULL) free (Ai) ;     \
    if (Ax != NULL) free (Ax) ;     \
    GrB_free (&A) ;                 \
    GrB_free (&C) ;                 \
}

static inline void func
(
    void *z,            // output, type uint8_t
    const void *aij,    // input, entry A(i,j), type uint8_t
    int64_t *indices,   // input, indices [i j]
    int n,              // 1 or 2
    const void *thunk   // inptut, a scalar, type uint8_t, unused
)
{
    uint8_t x = *((uint8_t *) aij) ;
    int64_t j = indices [1] ;           // row index ignored
    uint8_t result = x + (j & 0xFF) ;
    *((uint8_t *) z) = result ;
}

static inline void func2
(
    void *z,            // output, type uint8_t
    const void *aij,    // input, entry A(i,j), type uint8_t
    int64_t i,          // input
    int64_t j,          // input
    const void *thunk   // inptut, a scalar, type uint8_t, unused
)
{
    uint8_t x = *((uint8_t *) aij) ;
    uint8_t result = x + (j & 0xFF) ;
    *((uint8_t *) z) = result ;
}

static inline void func4
(
    void *z,            // output, type uint8_t
    const void *aij,    // input, entry A(i,j), type uint8_t
    int64_t i,          // input
    int64_t j,          // input
    const void *thunk   // inptut, a scalar, type uint8_t, unused
)
{
    uint8_t x = *((uint8_t *) aij) ;
    uint8_t result = x + (i & 0xFF) ;
    *((uint8_t *) z) = result ;
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;
    GrB_Matrix A = NULL, C = NULL ;
    LAGraph_Graph G = NULL ;
    uint8_t *Ax = NULL ;
    GrB_Index *Ap = NULL, *Ai = NULL ;
    size_t Ap_size = 0, Ai_size = 0, Ax_size = 0 ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    if (readproblem (&G, NULL,
        false, false, false, NULL, false, argc, argv) != 0) ERROR ;

    //--------------------------------------------------------------------------
    // convert to uint8
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols, nvals ;
    GrB_Matrix_nrows (&nrows, G->A) ;
    GrB_Matrix_ncols (&ncols, G->A) ;
    GrB_Matrix_nvals (&nvals, G->A) ;
    GrB_Matrix_new (&A, GrB_UINT8, nrows, ncols) ;
    GrB_Matrix_assign (A, NULL, NULL, G->A, GrB_ALL, nrows, GrB_ALL, ncols,
        NULL) ;
    LAGraph_Delete (&G, NULL) ;

    //--------------------------------------------------------------------------
    // export the A matrix as CSC
    //--------------------------------------------------------------------------

    char type_name [GxB_MAX_NAME_LEN] ;
    GrB_TRY (GxB_Matrix_type_name (type_name, A)) ;
    printf ("type: %s\n", type_name) ;
    GrB_TRY (GxB_Matrix_unpack_CSC (A, &Ap, &Ai, (void *) &Ax,
        &Ap_size, &Ai_size, &Ax_size, NULL, NULL, NULL)) ;
    GrB_free (&A) ;

    //--------------------------------------------------------------------------
    // apply using an idx-like operator with a single indices [2] array
    //--------------------------------------------------------------------------

    #define NTRIALS 10
    printf ("\nmethod 1:\n") ;

    double tic [2], tt ;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        uint8_t thunk = 0 ;     // unused
        LAGraph_Tic (tic, NULL) ;
        uint8_t *Cx = malloc (nvals * sizeof (uint8_t)) ;
        #pragma omp parallel for
        for (int64_t j = 0 ; j < ncols ; j++)
        {
            int64_t indices [2] ;
            indices [1] = j ;
            for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                indices [0] = Ai [p] ;      // read but not used
                // pretend to call a user-defined GrB_IndexUnaryOp
                func (&(Cx [p]), &(Ax [p]), indices, 2, &thunk) ;
            }
        }
        LAGraph_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
        printf ("time: %g\n", tt) ;
        free (Cx) ;
    }

    //--------------------------------------------------------------------------
    // apply using an idx-like operator with i and j separated
    //--------------------------------------------------------------------------

    printf ("\nmethod 2:\n") ;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        uint8_t thunk = 0 ;     // unused
        LAGraph_Tic (tic, NULL) ;
        uint8_t *Cx = malloc (nvals * sizeof (uint8_t)) ;
        #pragma omp parallel for
        for (int64_t j = 0 ; j < ncols ; j++)
        {
            for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // pretend to call a user-defined GrB_IndexUnaryOp
                func2 (&(Cx [p]), &(Ax [p]), Ai [p], j, &thunk) ;
            }
        }
        LAGraph_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
        printf ("time: %g\n", tt) ;
        free (Cx) ;
    }

    //--------------------------------------------------------------------------
    // apply using an idx-like operator with i and j separated, no load of i
    //--------------------------------------------------------------------------

    printf ("\nmethod 3:\n") ;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        uint8_t thunk = 0 ;     // unused
        LAGraph_Tic (tic, NULL) ;
        uint8_t *Cx = malloc (nvals * sizeof (uint8_t)) ;
        #pragma omp parallel for
        for (int64_t j = 0 ; j < ncols ; j++)
        {
            for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // pretend to call a user-defined GrB_IndexUnaryOp
                func2 (&(Cx [p]), &(Ax [p]), 0, j, &thunk) ;
            }
        }
        LAGraph_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
        printf ("time: %g\n", tt) ;
        free (Cx) ;
    }

    //--------------------------------------------------------------------------
    // apply using an idx-like operator with i and j separated, no load of i
    //--------------------------------------------------------------------------

    printf ("\nmethod 3b:\n") ;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        uint8_t thunk = 0 ;     // unused
        LAGraph_Tic (tic, NULL) ;
        uint8_t *Cx = malloc (nvals * sizeof (uint8_t)) ;
        #pragma omp parallel for
        for (int64_t j = 0 ; j < ncols ; j++)
        {
            for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // pretend to call a user-defined GrB_IndexUnaryOp
                // func2 (&(Cx [p]), &(Ax [p]), 0, j, &thunk) ;
                Cx [p] = Ax [p] + (j & 0xFF) ;
            }
        }
        LAGraph_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
        printf ("time: %g\n", tt) ;
        free (Cx) ;
    }

    //--------------------------------------------------------------------------
    // apply using an idx-like operator with i and j separated, require i and j
    //--------------------------------------------------------------------------

    printf ("\nmethod 4:\n") ;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        uint8_t thunk = 0 ;     // unused
        LAGraph_Tic (tic, NULL) ;
        uint8_t *Cx = malloc (nvals * sizeof (uint8_t)) ;
        #pragma omp parallel for
        for (int64_t j = 0 ; j < ncols ; j++)
        {
            for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // pretend to call a user-defined GrB_IndexUnaryOp
                func4 (&(Cx [p]), &(Ax [p]), Ai [p], j, &thunk) ;
            }
        }
        LAGraph_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
        printf ("time: %g\n", tt) ;
        free (Cx) ;
    }

    //--------------------------------------------------------------------------
    // apply using an idx-like operator with i and j separated, require i and j
    //--------------------------------------------------------------------------

    printf ("\nmethod 5:\n") ;
    for (int trial = 0 ; trial < NTRIALS ; trial++)
    {
        uint8_t thunk = 0 ;     // unused
        LAGraph_Tic (tic, NULL) ;
        uint8_t *Cx = malloc (nvals * sizeof (uint8_t)) ;
        #pragma omp parallel for
        for (int64_t j = 0 ; j < ncols ; j++)
        {
            for (int64_t p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                // pretend to call a user-defined GrB_IndexUnaryOp
                // func4 (&(Cx [p]), &(Ax [p]), Ai [p], j, &thunk) ;
                Cx [p] = Ax [p] + (Ai [p] & 0xFF) ;
            }
        }
        LAGraph_TRY (LAGraph_Toc (&tt, tic, NULL)) ;
        printf ("time: %g\n", tt) ;
        free (Cx) ;
    }

    LAGraph_FREE_ALL ;
    LAGraph_TRY (LAGraph_Finalize (msg)) ;
    return (0) ;
}

