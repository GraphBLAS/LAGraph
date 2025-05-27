//------------------------------------------------------------------------------
// LG_check_lcc: compute local clustering coefficients of a graph (simple)
//------------------------------------------------------------------------------

// LAGraph, (c) 2023 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Pascal Costanza, Intel, Belgium

//------------------------------------------------------------------------------

// An implementation of local clustering coefficient based on the descripion at
// https://en.wikipedia.org/wiki/Clustering_coefficient
// This method is for testing only, to check the result of other, faster methods.
// Do not benchmark this method; it is simple by design.

#define LG_FREE_ALL                         \
{                                           \
    GrB_free (&S) ;                         \
    GrB_free (&W) ;                         \
    GrB_free (&LCC) ;                       \
    LAGraph_Free ((void **) &Sp, msg) ;     \
    LAGraph_Free ((void **) &Si, msg) ;     \
    LAGraph_Free ((void **) &Tp, msg) ;     \
    LAGraph_Free ((void **) &Ti, msg) ;     \
    LAGraph_Free ((void **) &vb, msg) ;     \
    LAGraph_Free ((void **) &vx, msg) ;     \
}

#include <stdlib.h>
#include "LG_internal.h"
#include "LG_test.h"
#include "LG_test.h"
#include "LG_Xtest.h"

// assumes that the indices array is sorted
GrB_Index find(const GrB_Index* indices, GrB_Index n, GrB_Index index) {
    GrB_Index i = 0, j = n, h;
    while (i < j) {
        h = (i+j) >> 1;
        if (indices[h] < index) {
            i = h + 1;
        } else {
            j = h;
        }
    }
    return i;
}

// computes how many elements are in the intersection of both sets
// assumes that both arrays are sorted and don't have any duplicates
GrB_Index intersection_size(
        GrB_Index* x, GrB_Index nx,
        GrB_Index* y, GrB_Index ny
) {
    GrB_Index n = 0;
    while ((nx > 0) && (ny > 0)) {
        if (y[0] > x[0]) {
            GrB_Index* tmp = x; x = y; y = tmp;
            GrB_Index ntmp = nx; nx = ny; ny = ntmp;
        }
        GrB_Index index = find(y, ny, x[0]);
        if ((index < ny) && (y[index] == x[0])) {
            n++;
            y += index+1;
            ny -= index+1;
        } else {
            y += index;
            ny -= index;
        }
        x += 1;
        nx -= 1;
    }
    return n;
}

// This is a slow version.
int LG_check_lcc(
    // outputs:
    GrB_Vector *coefficients,     // the local clustering coefficients
    // inputs
    LAGraph_Graph G,        // input graph
    char *msg
)
{
#if LAGRAPH_SUITESPARSE

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    bool undirected = (G->kind == LAGraph_ADJACENCY_UNDIRECTED) ||
                      ((G->kind == LAGraph_ADJACENCY_DIRECTED) &&
                       (G->is_symmetric_structure == LAGraph_TRUE)) ;
    bool directed = !undirected ;

    GrB_Matrix S = NULL, T = NULL, W = NULL ;
    GrB_Index *Sp = NULL, *Si = NULL, *Tp = NULL, *Ti = NULL ;
    void *Sx = NULL , *Tx = NULL ;
    int8_t *vb = NULL ;
    double *vx = NULL ;
    GrB_Vector LCC = NULL ;
    GrB_Index n, ncols ;
    LG_ASSERT (coefficients != NULL, GrB_NULL_POINTER) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ; //set n to number of rows
    GRB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_ASSERT (n == ncols, GrB_INVALID_OBJECT) ;

    GRB_TRY (GrB_Vector_new (&LCC, GrB_FP64, n)) ;

    GrB_Matrix A = G->A ;

    GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n, n)) ;
    GRB_TRY(GrB_assign(S, A, GrB_NULL, true, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S)) ;
    if (G->nself_edges != 0) {
        GRB_TRY(GrB_select(S, GrB_NULL, GrB_NULL, GrB_OFFDIAG, S, 0, GrB_NULL)) ;
    }

    if (directed)
    {
        GRB_TRY(GrB_Matrix_new(&W, GrB_BOOL, n, n)) ;
        GRB_TRY(GrB_eWiseAdd(W, GrB_NULL, GrB_NULL, GrB_ONEB_BOOL, S, S, GrB_DESC_T1)) ;
        T = W ;
    } else {
        T = S ;
    }

    GrB_Index Sp_size, Si_size, Sx_size ;
    bool Siso ;
    GRB_TRY(GxB_Matrix_unpack_CSR(S, &Sp, &Si, &Sx, &Sp_size, &Si_size, &Sx_size, &Siso, NULL, GrB_NULL));
    LAGraph_Free (&Sx, msg) ;

    GrB_Index Tp_size, Ti_size, Tx_size ;
    bool Tiso ;
    if (directed) {
        GRB_TRY(GxB_Matrix_unpack_CSR(T, &Tp, &Ti, &Tx, &Tp_size, &Ti_size, &Tx_size, &Tiso, NULL, GrB_NULL));
        LAGraph_Free (&Tx, msg) ;
    } else {
        Tp = Sp; Tp_size = Sp_size;
        Ti = Si; Ti_size = Si_size;
    }

    LAGraph_Calloc ((void **) &vb, n, sizeof (int8_t), msg) ;
    LAGraph_Malloc ((void **) &vx, n, sizeof (double), msg) ;

    int64_t i;
    GrB_Index nvals = 0 ;

#pragma omp parallel for schedule(dynamic) reduction(+ : nvals)
    for (i = 0; i < n; i++) {
        GrB_Index *neighbors = Ti + Tp[i];
        GrB_Index k = Tp[i + 1] - Tp[i];
        if (k < 2) continue;

        GrB_Index j, esum = 0;
        for (j = 0; j < k; j++) {
            GrB_Index e = neighbors[j];
            GrB_Index *links = Si + Sp[e];
            GrB_Index l = Sp[e + 1] - Sp[e];
            if (l == 0) continue;
            if (undirected) {
                l = find(links, l, e);
            }
            esum += intersection_size(neighbors, k, links, l);
        }

        if (esum == 0) continue;

        if (undirected) {
            esum *= 2;
        }
        vb[i] = true;
        vx[i] = ((double) esum) / ((double) (k * (k - 1)));
        nvals++;
    }

    if (!directed)
    {
        // do not free these; they are aliases to Sp and Si:
        Tp = NULL ;
        Ti = NULL ;
    }

    GRB_TRY(GxB_Vector_pack_Bitmap(LCC, &vb, (void**)&vx, n*sizeof(int8_t), n*sizeof(double), false, nvals, GrB_NULL)) ;

    *coefficients = LCC ; LCC = NULL ;
    LG_FREE_ALL;
    return (GrB_SUCCESS);
#else
    return (GrB_NOT_IMPLEMENTED) ;
#endif
}
