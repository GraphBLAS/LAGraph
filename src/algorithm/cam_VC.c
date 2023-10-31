#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "LG_internal.h"

GrB_Info LG_VertexCover_cam(
    const LAGraph_Graph G,
    GrB_Vector *VC_out
)
{
    GrB_Info info;
    GrB_Index nrows, ncols;

    // Get adjacency matrix from the graph
    GrB_Matrix A;
    GrB_Matrix_dup(&A, G->A); // Duplicate the adjacency matrix to avoid modifying the original

    GrB_Matrix_nrows(&nrows, A);
    GrB_Matrix_ncols(&ncols, A);
    GrB_Index n = nrows;   // Number of vertices in G

    GrB_Vector VC;
    GrB_Vector_new(&VC, GrB_BOOL, n); // Initialize vertex cover vector to all false (or 0)
    GrB_assign(VC, NULL, NULL, 0, GrB_ALL, n, NULL);

    printf("Newly created VC vector\n");
    GxB_print(VC, GxB_COMPLETE);

    while (true)
    {
        // 1. Compute the degree of each vertex
        GrB_Vector degree;
        GrB_Vector_new(&degree, GrB_UINT64, n); // Initialize to 0

        GrB_Vector ones;
        GrB_Vector_new(&ones, GrB_UINT64, n);
        GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL);

        GrB_vxm(degree, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_UINT64, ones, A, GrB_DESC_T1); // Transpose A for matrix-vector multiplication
        GrB_Vector_free(&ones);

        // 2. Identify the vertex v with the maximum degree
        GrB_Index maxIndex;
        uint64_t maxValue = 0;
        for (GrB_Index i = 0; i < n; i++)
        {
            uint64_t value;
            GrB_Vector_extractElement(&value, degree, i);
            if (value > maxValue)
            {
                maxValue = value;
                maxIndex = i;
            }
        }

        // Check for termination: If the max degree is 0, then all edges are covered
        if (maxValue == 0)
        {
            GrB_Vector_free(&degree);
            break;
        }

        // 3. Update VC to include the vertex with maximum degree
        GrB_Vector_setElement(VC, (bool)true, maxIndex);


        // 4. Remove the edges connected to vertex v in the adjacency matrix
        for (GrB_Index j = 0; j < n; j++)
        {
            GrB_Matrix_setElement(A, (uint64_t)0, maxIndex, j); // Zero out the row
            GrB_Matrix_setElement(A, (uint64_t)0, j, maxIndex); // Zero out the column
        }

        GrB_Vector_free(&degree);
    }

    // VC now holds the vertices that are part of the approximate vertex cover.
    *VC_out = VC; // Pass the resulting vertex cover back to the caller

    GrB_Matrix_free(&A);

    return GrB_SUCCESS;
}
