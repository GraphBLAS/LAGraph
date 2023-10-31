
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#define LG_FREE_WORK        \
{                           \
    GrB_free (&frontier);   \
}

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&l_parent);   \
    GrB_free (&l_level);    \
}

#include "LG_internal.h"

GrB_Info LG_BreadthFirstSearch_cam
(
    GrB_Vector    *parents,
    const LAGraph_Graph G,
    GrB_Index      src
)
{
    GrB_Matrix A = G->A;

    GrB_Index N;
    GrB_Matrix_nrows(&N, A);

    GrB_Vector_new(parents, GrB_UINT64, N);
    GrB_Vector_setElement(*parents, src, src);

    GrB_Vector wavefront;
    GrB_Vector_new(&wavefront, GrB_UINT64, N);
    GrB_Vector_setElement(wavefront, 1UL, src);

    GrB_Index nvals;
    GrB_Vector_nvals(&nvals, wavefront);

    while (nvals > 0)
    {
        GrB_apply(wavefront, GrB_NULL, GrB_NULL, GrB_ROWINDEX_INT64, wavefront, 0UL, GrB_NULL);

        GrB_vxm(wavefront, *parents, GrB_NULL, GrB_MIN_FIRST_SEMIRING_UINT64, wavefront, A, GrB_DESC_RSC);

        GrB_apply(*parents, GrB_NULL, GrB_PLUS_UINT64, GrB_IDENTITY_INT64, wavefront, GrB_NULL);

        GrB_Vector_nvals(&nvals, wavefront);
    }

    GrB_free(&wavefront);

    // revert parents to 1-based indexing

    printf("Printing parent vector after algorithm\n");
    GxB_print(*parents, GxB_COMPLETE);

    return GrB_SUCCESS;

}