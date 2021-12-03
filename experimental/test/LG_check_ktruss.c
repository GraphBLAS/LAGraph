//------------------------------------------------------------------------------
// LG_check_ktruss: construct the ktruss of a graph (simple method)
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

// A very slow, bare-bones ktruss method.  This method is for testing only, to
// check the result of other, faster methods.  Do not benchmark this method; it
// is slow and simple by design.  G->A must be symmetric, with no entries on
// its diagonal.

#define LAGraph_FREE_WORK                       \
{                                               \
    LAGraph_Free ((void **) &Cp) ;              \
    LAGraph_Free ((void **) &Cj) ;              \
    LAGraph_Free ((void **) &Cx) ;              \
    LAGraph_Free ((void **) &Ax) ;              \
}

#define LAGraph_FREE_ALL                        \
{                                               \
    LAGraph_FREE_WORK ;                         \
    GrB_free (&C) ;                             \
}

#include "LG_internal.h"
#include "LG_test.h"

int LG_check_ktruss
(
    // output
    GrB_Matrix *C_handle,   // the ktruss of G->A, of type GrB_UINT32
    // input
    LAGraph_Graph G,        // the structure of G->A must be symmetric
    uint32_t k,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Matrix C = NULL ;
    GrB_Index *Cp = NULL, *Cj = NULL ;
    uint32_t *Cx = NULL ;
    void *Ax = NULL ;
    GrB_Index n, ncols, Cp_len, Cj_len, Cx_len, nvals1, nvals2 ;
    LG_CHECK (C_handle == NULL, -1003, "C_handle is NULL") ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1002, "graph is invalid") ;
    LG_CHECK (G->ndiag != 0, -104, "G->ndiag must be zero") ;
    if (G->kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
       (G->kind == LAGRAPH_ADJACENCY_DIRECTED &&
        G->A_structure_is_symmetric == LAGRAPH_TRUE))
    {
        // the structure of A is known to be symmetric
        ;
    }
    else
    {
        // A is not known to be symmetric
        LG_CHECK (false, -1005, "G->A must be symmetric") ;
    }
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_CHECK (n != ncols, -1001, "A must be square") ;

    //--------------------------------------------------------------------------
    // export G->A in CSR form and discard its values
    //--------------------------------------------------------------------------

    size_t typesize ;
    LAGraph_TRY (LG_check_export (G, &Cp, &Cj, &Ax, &Cp_len, &Cj_len, &Cx_len,
        &typesize, msg)) ;
    LAGraph_Free ((void **) &Ax) ;

    //--------------------------------------------------------------------------
    // allocate Cx
    //--------------------------------------------------------------------------

    Cx = (uint32_t *) LAGraph_Malloc (Cx_len, sizeof (uint32_t)) ;
    LG_CHECK (Cx == NULL, GrB_OUT_OF_MEMORY, "out of memory") ;

    //--------------------------------------------------------------------------
    // construct the k-truss of G->A
    //--------------------------------------------------------------------------

    while (true)
    {

        //----------------------------------------------------------------------
        // compute the # of triangles incident on each edge of C
        //----------------------------------------------------------------------

        // masked dot-product method: C{C}=C*C' using the PLUS_ONE semiring
        #if !defined ( COVERAGE )
        #pragma omp parallel for schedule(dynamic,1024)
        #endif
        for (int64_t i = 0 ; i < n ; i++)
        {
            // for each entry in C(i,:)
            for (int64_t p = Cp [i] ; p < Cp [i+1] ; p++)
            {
                const int64_t j = Cj [p] ;
                uint32_t cij = 0 ;
                // cij += C(i,:) * C(j,:)'
                int64_t p1 = Cp [i] ;
                int64_t p1_end = Cp [i+1] ;
                int64_t p2 = Cp [j] ;
                int64_t p2_end = Cp [j+1] ;
                while (p1 < p1_end && p2 < p2_end)
                {
                    int64_t j1 = Cj [p1] ;
                    int64_t j2 = Cj [p2] ;
                    if (j1 < j2)
                    { 
                        // C(i,j1) appears before C(j,j2)
                        p1++ ;
                    }
                    else if (j2 < j1)
                    { 
                        // C(j,j2) appears before C(i,j1)
                        p2++ ;
                    }
                    else // j1 == j2
                    { 
                        // C(i,j1) and C(j,j1) are the next entries to merge
                        cij++ ;
                        p1++ ;
                        p2++ ;
                    }
                }
                Cx [p] = cij ;
            }
        }

        //----------------------------------------------------------------------
        // import C in CSR form
        //----------------------------------------------------------------------

        GrB_TRY (GrB_Matrix_import_UINT32 (&C, GrB_UINT32, n, n,
            Cp, Cj, Cx, Cp_len, Cj_len, Cx_len, GrB_CSR_FORMAT)) ;
        GrB_TRY (GrB_Matrix_nvals (&nvals1, C)) ;

        //----------------------------------------------------------------------
        // keep entries >= k-2 and check for convergence
        //----------------------------------------------------------------------

        GrB_TRY (GrB_select (C, NULL, NULL, GrB_VALUEGE_UINT32, C, k-2, NULL)) ;
        GrB_TRY (GrB_Matrix_nvals (&nvals2, C)) ;
        if (nvals1 == nvals2)
        {
            // C is now the k-truss of G->A
            LAGraph_FREE_WORK ;
            (*C_handle) = C ;
            return (GrB_SUCCESS) ;
        }

        //----------------------------------------------------------------------
        // export C in CSR form for the next iteration and free it
        //----------------------------------------------------------------------

        GrB_TRY (GrB_Matrix_export_UINT32 (Cp, Cj, Cx,
            &Cp_len, &Cj_len, &Cx_len, GrB_CSR_FORMAT, C)) ;
        GrB_TRY (GrB_free (&C)) ;
    }
}

