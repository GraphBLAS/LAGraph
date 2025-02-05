#include "LG_internal.h"
#include "LAGraphX.h"
// add this algorithm to LAGraphX.h

int LAGraph_independent_set
(
    // output
    GrB_Vector *C,
    // add number of colors as output

    // input
    LAGraph_Graph G,
    char *msg
)
{
    bool verbose = false;

    GrB_Index n;
    GrB_Matrix_nrows(&n, G->A);

    /* initialize C to vector of 0s */
    // turn this into a local variable
    GRB_TRY(GrB_Vector_new(C, GrB_INT32, n));
    GrB_assign(*C, GrB_NULL, GrB_NULL, 0, GrB_ALL, n, GrB_NULL);

    // lg_set_format_hint -> bitmap

    /* weights initialized randomly
    *  seed of 20 was chosen arbitrarily */   
    GrB_Vector weight = NULL;
    GRB_TRY(GrB_Vector_new(&weight, GrB_UINT64, n));
    GRB_TRY(GrB_assign (weight, NULL, NULL, 0, GrB_ALL, n, NULL));
    LAGraph_Random_Seed(weight, 20, msg);
    /* DEBUG PRINT */ if (verbose) { printf("\n[ DEBUG ] weight vector\n"); GxB_print(weight, 3); }

    GrB_Vector in_curr_subset = NULL;
    GRB_TRY(GrB_Vector_new(&in_curr_subset, GrB_INT32, n)); // bool

    GrB_Vector max_weights = NULL;
    GRB_TRY(GrB_Vector_new(&max_weights, GrB_UINT64, n));

    GrB_Scalar reduced_scalar = NULL;
    GrB_Scalar_new(&reduced_scalar, GrB_INT32);

    GrB_Scalar color_scalar = NULL;
    GrB_Scalar_new(&color_scalar, GrB_INT32);
    
    /* algorithm start */
    for (int color = 1; color < n+1; color++) {
        /* vxm - find maximum of all neighboring weights */
        GRB_TRY(GrB_vxm(max_weights, GrB_NULL, GrB_NULL, GrB_MAX_TIMES_SEMIRING_UINT64, weight, G->A, GrB_NULL));
        /* DEBUG PRINT */ if (verbose) { printf("ITERATION %d\n[DEBUG] max_weights\n", color); GxB_print(max_weights, 3); }

        /* eWiseAdd - 1 if current weight > max neighboring weight */
        GRB_TRY(GrB_eWiseAdd(in_curr_subset, GrB_NULL, GrB_NULL, GrB_GT_UINT64, weight, max_weights, GrB_NULL));
        /* DEBUG PRINT */ if (verbose) { printf("[DEBUG] in_curr_subset\n"); GxB_print(in_curr_subset, 3); }

        /* reduce - add up all entries in in_curr_subset - if 0, break */
        GRB_TRY(GrB_reduce(reduced_scalar, GrB_NULL, GrB_PLUS_INT32, in_curr_subset, GrB_NULL));
        int32_t reduced_scalar_int;
        GRB_TRY(GrB_Scalar_extractElement(&reduced_scalar_int, reduced_scalar));
        /* DEBUG PRINT */ if (verbose) { printf("[DEBUG] reduced_scalar_int: %d\n", reduced_scalar_int); }
        if (reduced_scalar_int == 0) { break; }

        /* assign - write current color to C vector according to in_curr_subset mask */
        GrB_Scalar_setElement(color_scalar, color);
        GrB_assign(*C, in_curr_subset, GrB_NULL, color_scalar, GrB_ALL, n, GrB_NULL);
        if(verbose) { printf("[DEBUG] C, current color: %d\n", color); GxB_print(*C, 3);}
        
        /* assign - write 0 to weight according to in_curr_subset mask */
        GrB_assign(weight, in_curr_subset, GrB_NULL, 0, GrB_ALL, n, GrB_NULL);
        /* DEBUG PRINT */ if (verbose) { printf("[DEBUG] new weight\n"); GxB_print(weight, 3);}
    }

    return (GrB_SUCCESS) ;
}
