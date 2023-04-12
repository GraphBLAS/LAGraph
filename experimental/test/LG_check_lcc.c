//------------------------------------------------------------------------------
// LG_check_lcc: compute local clustering coefficients of a graph (simple)
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
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

#define LG_FREE_ALL                             \
{                                               \
    LAGraph_Free ((void **) &column, msg) ;     \
    free(neighbors) ;                           \
    free(links) ;                               \
}

#include <stdlib.h>
#include "LG_internal.h"
#include "LG_test.h"
#include "LG_test.h"
#include "LG_Xtest.h"

int comp_GrB_Index(const void *x, const void *y) {
    GrB_Index i = *(const GrB_Index*)x;
    GrB_Index j = *(const GrB_Index*)y;
    if (i < j) return -1;
    if (i > j) return +1;
    return 0;
}

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

void exclude(GrB_Index* indices, GrB_Index *n, GrB_Index element) {
    GrB_Index index = find(indices, *n, element);
    if ((index < *n) && (indices[index] == element)) {
        *n = *n - 1;
        memcpy(&indices[index], &indices[index+1], (*n-index)*sizeof(GrB_Index));
    }
}

GrB_Index intersection_size(
        GrB_Index* x, GrB_Index nx,
        GrB_Index* y, GrB_Index ny
) {
    GrB_Index n = 0;
    while ((nx > 0) && (ny > 0)) {
        if (y[0] < x[0]) {
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

int LG_check_lcc(
    // outputs:
    GrB_Vector *coefficients,     // the local clustering coefficients
    // inputs
    LAGraph_Graph G,        // input graph
    char *msg
) {

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Vector column = NULL ;
    GrB_Index *neighbors = NULL, *links = NULL ;
    GrB_Index n, ncols ;
    LG_ASSERT (coefficients != NULL, GrB_NULL_POINTER) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ; //set n to number of rows
    GRB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_ASSERT (n == ncols, GrB_INVALID_OBJECT) ;

    GRB_TRY (GrB_Vector_new (coefficients, GrB_FP64, n)) ;

    GRB_TRY (GrB_Vector_new (&column, GrB_BOOL, n)) ;

    GrB_Index e, i, j, k, nn, nl, esum ;
    neighbors = (GrB_Index*)malloc(n*sizeof(GrB_Index)) ;
    links = (GrB_Index*)malloc(n*sizeof(GrB_Index)) ;
    LG_ASSERT (neighbors != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (links != NULL, GrB_NULL_POINTER) ;

    GrB_Matrix A = G->A;

    if ((G->kind == LAGraph_ADJACENCY_UNDIRECTED) ||
            ((G->kind == LAGraph_ADJACENCY_UNDIRECTED) &&
            G->is_symmetric_structure == LAGraph_TRUE)) {

        for (i = 0; i < n; i++) {

            GRB_TRY (GrB_extract (column, GrB_NULL, GrB_NULL, A, GrB_ALL, n, i, GrB_DESC_R)) ;
            nn = n;
            GRB_TRY (GrB_Vector_extractTuples (neighbors, (bool*)NULL, &nn, column)) ;
            if (nn == 0) continue;
            qsort(neighbors, nn, sizeof(GrB_Index), comp_GrB_Index);
            exclude(neighbors, &nn, i);
            k = nn;
            if (k < 2) continue;

            esum = 0;
            for (j = 0; j < nn; j++) {
                e = neighbors[j];
                GRB_TRY (GrB_extract (column, GrB_NULL, GrB_NULL, A, GrB_ALL, n, e, GrB_DESC_R));
                nl = n;
                GRB_TRY (GrB_Vector_extractTuples (links, (bool*)NULL, &nl, column));
                if (nl == 0) continue;
                qsort(links, nl, sizeof(GrB_Index), comp_GrB_Index);
                nl = find(links, nl, e);
                esum += intersection_size(neighbors, nn, links, nl);
            }
            GRB_TRY (GrB_Vector_setElement (*coefficients, ((double)(2*esum))/((double)(k*(k-1))), i));
        }
    } else {

        for (i = 0; i < n; i++) {
            GRB_TRY (GrB_extract (column, GrB_NULL, GrB_NULL, A, GrB_ALL, n, i, GrB_DESC_R)) ;
            GRB_TRY (GrB_extract (column, GrB_NULL, GxB_ANY_UINT32, A, GrB_ALL, n, i, GrB_DESC_T0)) ;
            nn = n;
            GRB_TRY (GrB_Vector_extractTuples (neighbors, (bool*)NULL, &nn, column)) ;
            if (nn == 0) continue;
            qsort(neighbors, nn, sizeof(GrB_Index), comp_GrB_Index);
            exclude(neighbors, &nn, i);
            k = nn;
            if (k < 2) {
                continue;
            }

            esum = 0;
            for (j = 0; j < nn; j++) {
                e = neighbors[j];
                GRB_TRY (GrB_extract (column, GrB_NULL, GrB_NULL, A, GrB_ALL, n, e, GrB_DESC_R));
                nl = n;
                GRB_TRY (GrB_Vector_extractTuples (links, (bool*)NULL, &nl, column));
                if (nl == 0) continue;
                qsort(links, nl, sizeof(GrB_Index), comp_GrB_Index);
                exclude(links, &nl, e);
                esum += intersection_size(neighbors, nn, links, nl);
            }
            GRB_TRY (GrB_Vector_setElement (*coefficients, ((double)esum)/((double)(k*(k-1))), i));
        }
    }

    LG_FREE_ALL;
    GRB_TRY (GrB_wait (*coefficients, GrB_MATERIALIZE));
    return (GrB_SUCCESS);
}
