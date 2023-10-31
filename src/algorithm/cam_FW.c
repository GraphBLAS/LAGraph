#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

// #define LG_FREE_ALL        \
// {                          \
//     GrB_free(D) ;         \     \    
// }


#include "LG_internal.h"

GrB_Info LG_FloydWarshall_cam 
(
    const LAGraph_Graph G,
    GrB_Matrix *D_out
) 
{
    GrB_Info info;
    char *msg = NULL;
    GrB_Index nrows, ncols;

    GrB_Matrix D = NULL;
    GrB_Matrix D_col_k, D_row_k, D_k = NULL;

    GrB_Matrix A = G->A;

    GrB_Matrix_nrows(&nrows, A);
    GrB_Matrix_ncols (&ncols, A);
    LG_ASSERT_MSG (nrows == ncols, -1002, "A must be square") ;
    GrB_Index n = nrows;   // n = size of N x N matrix A    
    
    // Initialize D_0 to A
    GrB_Matrix_dup(&D, A);

    // Create temporary vectors for the kth column and row of D
    GrB_Matrix_new(&D_col_k, GrB_FP64, n, (GrB_Index) 1);
    GrB_Matrix_new(&D_row_k, GrB_FP64, (GrB_Index) 1, n);

    for (GrB_Index k = 0; k < n; k++)
    {
        printf("Iteration %i\n", k);
        // GxB_print (D, GxB_COMPLETE);
        // Extract kth column and row of D
        GrB_extract(D_col_k, NULL, NULL, D, GrB_ALL, n, &k, 1, NULL);
        GrB_extract(D_row_k, NULL, NULL, D, &k, 1, GrB_ALL, n, NULL);
        GxB_print (D_row_k, GxB_COMPLETE);
        GxB_print (D_col_k, GxB_COMPLETE);

        // GrB_eWiseMult(D_k, NULL, NULL, GrB_MIN_PLUS_SEMIRING_FP64, D_col_k, D_row_k, NULL);
        GrB_mxm(D, NULL, GrB_MIN_FP64, GrB_MIN_PLUS_SEMIRING_FP64, D_col_k, D_row_k, NULL);
        GxB_print (D, GxB_COMPLETE);

        // D_k = D_{k - 1} .min ( D_{k - 1}(:, k) min.+ D_{k - 1}(k, :) )
        // GrB_eWiseAdd(D, NULL, NULL, GrB_MIN_MONOID_FP64, D, D_k, NULL);
    }

    *D_out = D;

    // Cleanup
    GrB_free(&D_col_k);
    GrB_free(&D_row_k);
    GrB_free(&D_k);

    return GrB_SUCCESS;
}

